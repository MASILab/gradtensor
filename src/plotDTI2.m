function plotDTI2(ax,D,delta)
if nargin==1
    delta=1;
end
sz=size(D);
if length(sz)==2
    nx=1;ny=1;
elseif length(sz)==3
    nx=sz(3);ny=1;
elseif length(sz)==4
    nx=sz(3);ny=sz(4);
elseif length(sz)==5
    nx=sz(4);ny=sz(5);
end
disp(nx)
disp(ny)
n=size(D,3);
for i=1:nx
    for j=1:ny
        D=squeeze(D);
        [v,l]=eig(D(:,:,i,j));
        [X,Y,Z]=ellipsoid(0,0,0,l(1,1),l(2,2),l(3,3),10);
        sz=size(X);
        for x=1:sz(1)
            for y=1:sz(2)
                A=[X(x,y) Y(x,y) Z(x,y)]';
                A=v*A;
                X(x,y)=A(1);Y(x,y)=A(2);Z(x,y)=A(3);
            end
        end
        X=X+(i-1)*delta*2;
        Y=Y+(j-1)*delta*2;
        surf(ax,X,Y,Z);
        %xlim([-1 .001]);
        %ylim([-1 .001]);
        %zlim([-1 .001]);
        %set(h, 'FaceAlpha', 0.5)
        %alpha 0.5
        %if i==1 & j==1
         %   hold on
        end
    end
end

% axis equal;
% alpha 0.5;
% disp(xlim);
% disp(ylim);
% disp(zlim);
% xlim([-0.002 0.002]);
% ylim([-0.002 0.002]);
% zlim([1.0e-03 * -0.8   1.0e-03 * 0.8]);
% %view([90 0]);
% set(gca,'GridLine','none');
% %set(gca,'xlimMode',[-1 0.003]);
% set(gca,'XTick');
% set(gca,'YTick');
% set(gca,'ZTick');
% colormap(ax,[0.8 0.8 0.8])
% hold on
