
clear all

load('field7R0','PFF','omega','mesh','dim','index')



PFF1 = reshape(PFF(index.block1(:)),dim.dimY1,dim.dimX1);
PFF2 = reshape(PFF(index.block2(:)),dim.dimY2,dim.dimX2);
PFF3 = reshape(PFF(index.block3(:)),dim.dimY3,dim.dimX3);

AbsPFF1=abs(PFF1); AbsPFF2=abs(PFF2); AbsPFF3=abs(PFF3);

ArgPFF1=angle(PFF1); ArgPFF2=angle(PFF2); ArgPFF3=angle(PFF3);

%%


%
figure(1)
hold off

surf(mesh.X1,mesh.Y1,full(AbsPFF1));
hold on
surf(mesh.X2,mesh.Y2,full(AbsPFF2));
hold on
surf(mesh.X3,mesh.Y3,full(AbsPFF3));
hold on
surf(-mesh.X1,mesh.Y1,full(AbsPFF1));
hold on
surf(-mesh.X2,mesh.Y2,full(AbsPFF2));
hold on
surf(-mesh.X3,mesh.Y3,full(AbsPFF3));
view(0,90)
axis equal
maxAbsP=max([max(abs(AbsPFF1(:))),max(abs(AbsPFF2(:))),max(abs(AbsPFF3(:)))]);
minAbsP=min([min(abs(AbsPFF1(:))),min(abs(AbsPFF2(:))),min(abs(AbsPFF3(:)))]);
caxis([-maxAbsP maxAbsP]);
xlabel('[m]','FontSize',20,'FontWeight','bold','FontName','times')
ylabel('[m]','FontSize',20,'FontWeight','bold','FontName','times')
shading flat
colormap(jet(200));
colorbar
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',16,'FontWeight','bold')

figure(2)
hold off

surf(mesh.X1,mesh.Y1,full(ArgPFF1));
hold on
surf(mesh.X2,mesh.Y2,full(ArgPFF2));
hold on
surf(mesh.X3,mesh.Y3,full(ArgPFF3));
hold on
surf(-mesh.X1,mesh.Y1,full(ArgPFF1));
hold on
surf(-mesh.X2,mesh.Y2,full(ArgPFF2));
hold on
surf(-mesh.X3,mesh.Y3,full(ArgPFF3));
view(0,90)
axis equal
maxArg=max([max(abs(ArgPFF1(:))),max(abs(ArgPFF2(:))),max(abs(ArgPFF3(:)))]);
minArg=min([min(abs(ArgPFF1(:))),min(abs(ArgPFF2(:))),min(abs(ArgPFF3(:)))]);
caxis([-maxArg maxArg]);
xlabel('[m]','FontSize',20,'FontWeight','bold','FontName','times')
ylabel('[m]','FontSize',20,'FontWeight','bold','FontName','times')
shading flat
colormap(jet(200));
colorbar
set(gcf, 'color', [1 1 1])
set(gca,'FontSize',16,'FontWeight','bold')

%%
%

N=100; filename='mode1conf7R0.gif'

phase=[0:2*pi/(N-1):4*pi]

for ns=1:N
    
    
    PT1 = real(PFF1.*exp(1i*phase(ns)));
    PT2 = real(PFF2.*exp(1i*phase(ns)));
    PT3 = real(PFF3.*exp(1i*phase(ns)));
    
    figure(3)
    hold off
    
    surf(mesh.X1,mesh.Y1,full(PT1));
    hold on
    surf(mesh.X2,mesh.Y2,full(PT2));
    hold on
    surf(mesh.X3,mesh.Y3,full(PT3));
    hold on
    surf(-mesh.X1,mesh.Y1,full(PT1));
    hold on
    surf(-mesh.X2,mesh.Y2,full(PT2));
    hold on
    surf(-mesh.X3,mesh.Y3,full(PT3));
    view(0,90)
    %view(3)
    axis equal
    caxis([-maxAbsP maxAbsP]);
    %xlabel('[m]')
    %ylabel('[m]')
    shading flat
    colormap(jet(200));
    colorbar
     set(gcf, 'color', [1 1 1])
    
    %pause(1)
    
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,1280);
    
    if ns == 1;
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0);
    end
    
    
    
end

%}
