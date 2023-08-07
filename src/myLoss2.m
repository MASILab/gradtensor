function loss = myLoss2(Scorrupt,b,qii,L_mat2)

Lxx = L_mat2(1);
Lxy = L_mat2(2);
Lxz = L_mat2(3);

Lyx = L_mat2(4);
Lyy = L_mat2(5);
Lyz = L_mat2(6);


Lzx = L_mat2(7);
Lzy = L_mat2(8);
Lzz = L_mat2(9);

Lii2 = transpose([Lxx^2  Lxx*Lxy  Lxx*Lxz  Lxx*Lyx  Lxx*Lyy  Lxx*Lyz  Lxx*Lzx  Lxx*Lzy  Lxx*Lzz  Lxy^2 ...
    Lxy*Lxz  Lxy*Lyx  Lxy*Lyy  Lxy*Lyz  Lxy*Lzx  Lxy*Lzy  Lxy*Lzz  Lxz^2  Lxz*Lyx  Lxz*Lyy  ...
    Lxz*Lyz  Lxz*Lzx  Lxz*Lzy  Lxz*Lzz  Lyx^2  Lyx*Lyy  Lyx*Lyz  Lyx*Lzx  Lyx*Lzy  Lyx*Lzz  ...
    Lyy^2  Lyy*Lyz  Lyy*Lzx  Lyy*Lzy  Lyy*Lzz  Lyz^2  Lyz*Lzx  Lyz*Lzy  Lyz*Lzz  Lzx^2 ...
    Lzx*Lzy  Lzx*Lzz  Lzy^2  Lzy*Lzz  Lzz^2]);

% Lii2 = [L(1)^2  L(1)*L(2)  L(1)*L(3)  L(1)*L(4)  L(1)*L(5)  L(1)*L(6)  L(1)*L(7)  L(1)*L(8)  L(1)*L(9)  L(2)^2 ...
%     L(2)*L(3)  L(2)*L(4)  L(2)*L(5)  L(2)*L(6)  L(2)*L(7)  L(2)*L(8)  L(2)*L(9)  L(3)^2  L(3)*L(4)  L(3)*L(5)  ...
%     L(3)*L(6)  L(3)*L(7)  L(3)*L(8)  L(3)*L(9)  L(4)^2  L(4)*L(5)  L(4)*L(6)  L(4)*L(7)  L(4)*L(8)  L(4)*L(9)  ...
%     L(5)^2  L(5)*L(6)  L(5)*L(7)  L(5)*L(8)  L(5)*L(9)  L(6)^2  L(6)*L(7)  L(6)*L(8)  L(6)*L(9)  L(7)^2 ...
%     L(7)*L(8)  L(7)*L(9)  L(8)^2  L(8)*L(9)  L(9)^2]';
loss = norm(Scorrupt'- exp(-b.*qii*Lii2));
% loss = norm(log(Scorrupt)'- (-b.*qii*Lii2));