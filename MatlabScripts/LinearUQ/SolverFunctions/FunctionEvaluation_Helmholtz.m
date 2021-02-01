
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT(3){

% evaluate the response function values at integration points :

clearvars -except n1 n2 n3 n4 q1 q2 q3 q4 num
for vari4=1:q4
    for vari3=1:q3
        for vari2=1:q2
            for vari1=1:q1


close all
%clearvars -except vari1 vari2 vari3 vari4 Output1 Output2 n1 n2 n3 n4

prvi=n1;
drugi=n2;
tretji=n3;%[1.35 1.5]
cetrti=n4;%[0.9*4.73e-3 4.73e-3 ]

%%%%  Get Matrix A %%%%
configs=[8];  slope=0.0;

FlameType = 'NTAU' % 'NTAU', 'FTF'
FLAME = 'B'; % 'A', 'B'

shape=1; thin=1; 

MagRn=prvi(vari1); ArgRn=drugi(vari2);

fgain=tretji(vari3);  %Amp=0.1;


if strcmp(FLAME,'A'); tdelay=6.0e-3; end
if strcmp(FLAME,'B'); tdelay=cetrti(vari4); end



%file=['./RESULT_','_FlameType_',FlameType, ...
   % '_gain_',num2str(fgain),'_tdelay',num2str(tdelay),'_C',num2str(configs),'.mat'];



%=====

Helmholtz_Body

% fref (Hz)  and gref (Rad/s)

%save(file,'configs','fref','gref');

Output1(vari1,vari2,vari3,vari4)=real(omega)/2/pi;% angular frequency wr
Output2(vari1,vari2,vari3,vari4)=imag(omega);% growth rate wi

            end
        end
    end
end

 

% give a name to your function values data and save it in text file

fvalues=reshape(Output1,[num,1]);

fileID = fopen('FunctionValuesUNI.txt','wt');
fprintf(fileID,'%14.10f \n',fvalues');
fclose(fileID); 

% reading fvalues from your text file:
%[fvalues] = textread('FunctionValuesUNI.txt','%f');


%give a name to your function values data and save it 
save Outputs_config8_3333.mat Output1 Output2
disp('data successfully saved')



% } USER INPUT(3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





