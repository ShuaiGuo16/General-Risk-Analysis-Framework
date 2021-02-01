function [fref,gref] = Helmholtz_rf(solver,fit,MagRs,ArgRs,MagRn,ArgRn,alpha_damp,configs,s_init)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helmholtz solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%       solver:              Iteration algorithm ('Fixpoint' or 'Secant')
%       fit:                    rational function object     
%       MagRs:            magnitude of reflection coefficient at inlet
%       ArgRs:             phase of reflection coefficient at inlet
%       MagRn:            magnitude of reflection coefficient at outlet
%       ArgRn:             phase of reflection coefficient at outlet
%       alpha_damp:     system acousticdamping coefficient
%       configs:             configuration index
%       s_init:                initial value for iteration algorithm
% Outputs:
%       fref:                   modal frequency
%       gref:                  modal growth rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by C. F. Silva (TUM), Sept. 2016
% Email: camilo.f.silva.g@gmail.com
% Version: MATLAB R2015b
% Ref: [1] 
%      [2] Gustavsen.B and A.Semlyen, ��Rational approximation of 
%          frequency domain responses by vector fitting,�� IEEE 
%          Trans. Power Delivery, Vol. 14, No. 3, pp. 1052�C1061, 1999.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amp=0.1; slope=0.0; shape=1; thin=1; iterlim=400;

for cn=1:length(configs)
    
    conf=configs(cn);
    InitFVM
    [mesh, index] = setUpMesh(l,dim);
    
    if shape==0; cfield; else; cfield2; end
    
    % bound = boundary;
    
    [A,B,Axn,Axs] = getLinSyst(mesh, index, bound, dim);
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize w an eigenmode for faster calculations %%%
    s = s_init;
    
    %     fprintf('\n Initialize ready. Natural Freqency: %d \n\n',imag(s)/(2*pi));
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % Set values for FDF FLAME A
    flame.uref = 2.67; % m/s
    flame.Qref = 1940; % W
    flame.rhoref = 1.177; % kg/m3
    flame.gamma = 1.4; % [-]
    flame.alpha = 0.5; % Relaxation Parameter
    
    Reflexn = MagRn*exp(1i*ArgRn);
    Reflexs = MagRs*exp(1i*ArgRs);
    
    Zn=(1+Reflexn)/(1-Reflexn);
    Zs=(1+Reflexs)/(1-Reflexs);
    
    alphan = 1/(Zn*c_sound2);
    alphas = 1/(Zs*c_sound1);
    I=speye(length(A));
    switch solver
        case 'Fixpoint'
            lambda = s^2;
            iter=1; % Counter for relaxation
            ds=1;
            
            
            while abs(ds) > 1e-4 && iter<iterlim
                tic
                
                %omega=w;
                
                %%%%%%% Construct Flame Model %%%%%%%%
                
                flame.FTF = sum(fit.C./(s-fit.A));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Get nonlinear part H of matrix A due to flame
                H = sflame(mesh, index, dim, l, flame, thin);
                
                % Merge space discretization scheme and flame
                
                Anl = A + s*alphan*Axn + s*alphas*Axs - alpha_damp*s*I - H;
                
                
                % Solve the eigenvalue problem
                numEigValue = 1;% Number of modes to be calculated around zero, don't forget
                [PFF,lambda] = eigs(Anl,numEigValue,lambda);
                
                % Get first frequency from omega
                
                snew = sqrt(lambda);
                
                if imag(snew) < 0; snew = -sqrt(lambda); end  % it is important to chose physical eigenvalue
                
                ds = imag(snew) - imag(s);
                s = flame.alpha*snew + (1-flame.alpha)*s;
                lambda = s^2;
                
                %fprintf('Iteration %d: w = %d, time: %d\n',iter,w,toc);
                
                iter=iter+1;
                
            end
        case 'Secant'
            s1=s;
            s2=0.9*s;
            svec =[s1;s2];
            lambdavec = svec.^2;
            iter=1; % Counter for relaxation
            ds=1;
            while 1 %abs(ds) > 1e-2 && iter < iterlim
                tic
                
                %%%%%%% Construct Flame Model %%%%%%%%%
                
                flame.FTF = sum(fit.C./(svec(1)-fit.A));
                H1 = sflame(mesh, index, dim, l, flame, thin);
                flame.FTF = sum(fit.C./(svec(2)-fit.A));
                H2 = sflame(mesh, index, dim, l, flame, thin);

                if iter ==1
                    
                    
                    % Merge space discretization scheme and flame
                    
                    Anl1 = A + svec(1)*alphan*Axn + svec(1)*alphas*Axs - alpha_damp*svec(1)*I - H1;
                    Anl2 = A + svec(2)*alphan*Axn + svec(2)*alphas*Axs - alpha_damp*svec(2)*I - H2;
                    
                    L1 = lambdavec(1)*I-Anl1;
                    L2 = lambdavec(2)*I-Anl2;
                    
                    
                    YY1=eigs(L1,1,0);
                    YY2=eigs(L2,1,0);
                else
                    YY2=YY1;
                    Anl1 = A + svec(1)*alphan*Axn + svec(1)*alphas*Axs - alpha_damp*svec(1)*I - H1;
                    L1 = lambdavec(1)*I-Anl1;
                    YY1=eigs(L1,1,0);
                end
                
                lambdanew = lambdavec(1) - YY1*(lambdavec(1)-lambdavec(2))/(YY1-YY2);
                
                % Get first frequency from omega
                
                snew = sqrt(lambdanew);
                
                if imag(snew) < 0; snew = -sqrt(lambdanew); end  % it is important to chose physical eigenvalue
                
                ds0=ds;
                
                ds = imag(snew) - imag(svec(1));
                
                if abs(ds) > 1e-4 && iter < iterlim
                    
                    if iter>10 &&   abs(ds)-abs(ds0) > 10
                        flame.alpha = flame.alpha/2;
                    end
                    svec(2) = svec(1);
                    svec(1) = flame.alpha*snew + (1-flame.alpha)*svec(1);
                    %                             w = s/1i;
                    
                    lambdavec = svec.^2;
%                     if rem(iter,10)==1
%                         fprintf('Iteration %d: f = %f%+fi, ds= %d\n',iter,imag(svec(1))/2/pi,real(svec(1))/2/pi,abs(ds));
%                     end
                    iter=iter+1;
                else
                    
                    if iter >= iterlim
                        s = 0; lambda = 0;
                    else
                        [PFF,YYerr]=eigs(L1,1,0);
                        s=svec(1);
                        lambda=s^2;
                    end
                    break
                    
                end
            end
    end
    omega = sqrt(-lambda);
    
    fref(cn)       = imag(s)/2/pi;
    gref(cn)       = real(s);
    
    
end

end