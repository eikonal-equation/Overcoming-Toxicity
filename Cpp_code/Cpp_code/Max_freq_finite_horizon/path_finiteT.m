function [xlist,ylist,policy_list]= path_finiteT(xc,yc,n,dt,Tf,A,choice,gamma,epsilon,r_ks)
% N = length(U);
h = 1/n;
xx = 0:h:1;

% velocity function for different r_k and r_s
fx = @(x,y,a)  x.*(1-x).*((1-y).*(r_ks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);

fy = @(x,y,a)  y.*(1-y).*(1 + (r_ks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);



tt = 0:dt:Tf; %number of time steps

tnum = length(tt);


%% Euler
ylist=zeros(tnum,1);
xlist=zeros(tnum,1);
ylist(1) = yc;
xlist(1) = xc;
a1_set=[];
a0_set=[];
policy_list=[];

for nn = 1:tnum-1
    
    yval = ylist(nn);
    xval = xlist(nn);

    kx = find(xlist(nn) <= xx, 1);
    ky = find(ylist(nn) <= xx, 1);

    A_bottom = A(:,:,nn);

    a1=A_bottom(ky-1,kx-1);a2=A_bottom(ky-1,kx);
    a3=A_bottom(ky,kx-1);a4=A_bottom(ky,kx);
    a_square = [a1,a2,a3,a4];
    switch choice
        case 'conservative'
            %super-conservative
            if sum(a_square)== 4
                policy = 1;

                a1_set=[a1_set,nn];

            else
                policy = 0;

                a0_set=[a0_set,nn];

            end
            policy_list=[policy_list,policy];
        case 'aggressive'
            %super-aggressive
            if sum(a_square) > 0
                policy = 1;

                a1_set=[a1_set,nn];

            else
                policy = 0;

                a0_set=[a0_set,nn];

            end
            policy_list=[policy_list,policy];
        case 'majority'
            %majority
            if sum(a_square)>=3
                policy = 1;

                a1_set=[a1_set,nn];

            else
                policy = 0;

                a0_set=[a0_set,nn];

            end
    end
    policy_list=[policy_list,policy];


    dxx = fx(xval,yval,policy);
    dyy = fy(xval,yval,policy);


    ylist(nn+1) = yval+dyy*dt;
    xlist(nn+1) = xval+dxx*dt;





end



end


