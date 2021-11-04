function plotDTI(D,delta)%, color)
%-fanDTasia ToolBox------------------------------------------------------------------
% This Matlab script is part of the fanDTasia ToolBox: a Matlab library for Diffusion 
% Weighted MRI (DW-MRI) Processing, Diffusion Tensor (DTI) Estimation, High-order 
% Diffusion Tensor Analysis, Tensor ODF estimation, Visualization and more.
%
% A Matlab Tutorial on DW-MRI can be found in:
% http://www.cise.ufl.edu/~abarmpou/lab/fanDTasia/tutorial.php
%
%-CITATION---------------------------------------------------------------------------
% If you use this software please cite the following work:
% A. Barmpoutis, B. C. Vemuri, T. M. Shepherd, and J. R. Forder "Tensor splines for 
% interpolation and approximation of DT-MRI with applications to segmentation of 
% isolated rat hippocampi", IEEE TMI: Transactions on Medical Imaging, Vol. 26(11), 
% pp. 1537-1546 
%
%-DESCRIPTION------------------------------------------------------------------------
% This function plots a 2D field of 3D tensors as ellipsoidal glyphs. The 3D tensors 
% must be in the form of 3x3 symmetric positive definite matrices. The field can 
% contain either a single tensor, or a row of tensors or a 2D field of tensors.
%
% NOTE: This function plots only 2nd-order tensors (i.e. traditional DTI). For higher
% order tensors (such as 4th-order tensor visualization) please use the plotTensors.m
%
%-USE--------------------------------------------------------------------------------
% example 1: plotDTI(D)
% where D is of size 3x3 or 3x3xN or 3x3xNxM
%
% example 2: plotDTI(D,delta)
% where delta is a scalar that controls the size 
% of a voxel in the field. Default: delta=1
%
%-DISCLAIMER-------------------------------------------------------------------------
% You can use this source code for non commercial research and educational purposes 
% only without licensing fees and is provided without guarantee or warrantee expressed
% or implied. You cannot repost this file without prior written permission from the 
% authors. If you use this software please cite the following work:
% A. Barmpoutis, B. C. Vemuri, T. M. Shepherd, and J. R. Forder "Tensor splines for 
% interpolation and approximation of DT-MRI with applications to segmentation of 
% isolated rat hippocampi", IEEE TMI: Transactions on Medical Imaging, Vol. 26(11), 
% pp. 1537-1546 
%
%-AUTHOR-----------------------------------------------------------------------------
% Angelos Barmpoutis, PhD
% Computer and Information Science and Engineering Department
% University of Florida, Gainesville, FL 32611, USA
% abarmpou at cise dot ufl dot edu
%------------------------------------------------------------------------------------
%ax1 = axes;
%ax2 = axes;
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
        surf(X,Y,Z)%, 'FaceColor', color , 'FaceAlpha', 0.3)%, 'EdgeColor', 'none');
        %hold on 
        %plot3(X,Y,Z)
        %surf(X,Y,Z, 'FaceAlpha', 0.5);
        %xlim([-1 .001]);
        %ylim([-1 .001]);
        %zlim([-1 .001]);
        %set(h, 'FaceAlpha', 0.5)
        %alpha 0.5
        if i==1 & j==1
            hold on
        end
    end
end
axis equal;

%xlim([1.0e-03 * -2   1.0e-03 * 2]);
%xlim([-0.0021 0.0021]);
%ylim([1.0e-03 * -1.5   1.0e-03 * 1.5]);
%zlim([1.0e-4 * -6.5   1.0e-04 * 6.5]);
%xlim([-0.013 0.013]);
%ylim([1.0e-03 * -0.5470 1.0e-03 * 0.5470]);
%zlim([1.0e-03 * -0.5202   1.0e-03 * 0.5202]);
%xlim([-7 7])
%ylim([-0.03 0.03])
%zlim([-0.06 0.06])
%view([-154 -4]);
set(gca,'GridLine','none');
%set(gca,'xlimMode',[-1 0.003]);
set(gca,'XTick');
set(gca,'YTick');
set(gca,'ZTick');
alpha 0.3
%shading interp
%colormap([0 0.57 0.57])
colormap(BWR2)
%lighting phong
%light
%light('Position',[0 0 1]);%,'Style','infinite','Color',[ 1.000 0.584 0.000]);
xlabel('X');
ylabel('Y');
zlabel('Z');
disp(xlim);
disp(ylim);
disp(zlim);
axis off
%hold on
%grid on
%fprintf(1,'\nIf you use plotDTI.m please cite the following work:\n');
%fprintf(1,'A. Barmpoutis, B. C. Vemuri, T. M. Shepherd, and J. R. Forder "Tensor splines for\n');
%fprintf(1,'interpolation and approximation of DT-MRI with applications to segmentation of\n');
%fprintf(1,'isolated rat hippocampi", IEEE TMI: Transactions on Medical Imaging, Vol. 26(11),\n');
%fprintf(1,'pp. 1537-1546\n');