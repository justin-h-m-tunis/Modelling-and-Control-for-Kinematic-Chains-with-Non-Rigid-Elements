%KDC final

%start params
    %n
        n = 2;
    %q (state variable)
        syms t
        q = sym('q',[2*n,6]) ;
        q_th = sym('qt',[n,1]);
        q_del = sym('qd',[n,6]);
        ind_th = mod(1:2*n,2)==1;
        ind_del = mod(1:2*n,2)==0;
        q(ind_th,1) = q_th;
        q(ind_th,2:6) = zeros(n,5);
        q(ind_del,:) = q_del;
        dq = sym('dq',[2*n,6]);
        dq_th = sym('dqt',[n,1]);
        dq_del = sym('dqd',[n,6]);
        dq(ind_th,1) = dq_th;
        dq(ind_del,:) = dq_del;
        dq(ind_th,2:6) = zeros(n,5);
        d2q = sym('d2q',[2*n,6]);
        d2q(ind_th,2:6) = zeros(n,5);
    %w_i
        w_i = sym('w',[n,3]);
        w_i(:,1:2) = zeros(n,2);
        w_i(:,3) = ones(n,1);
    %q_i (axis and tool location)
        l = sym('l',[n,1]);
        lengths = [sym(0);l];
        q_i = sym('qi',[n+1,3]);
        for i = 1:n+1
            q_i(i,1) = sum(lengths(1:i));
            q_i(i,2:3) = [0,0];
        end
    %m_th,I_th
        m_i = sym('m',[2*n,1]);
        I_i = sym('I',[3,3,2*n]);
    %E,IA_xyr
        syms E G A IA_y IA_z IA_r tan_del;
    %twist order (order to apply the bending to beam)
        twist_order = [2,3,5,6,4]; % XXaxialXX, bend_y, bend_z, torsion_y, torsion_z, twist_x
    %g_0
        r = sym('r',[3,n]);
        g_0 = sym('g0',[4 4 2*n]);
        for i = 1:n
            g_0(:,:,2*i-1) = [eye(3),q_i(i,:).'+r(:,i); 0 0 0 1];
            g_0(:,:,2*i) = [eye(3),q_i(i+1,:).'; 0 0 0 1];
        end
%derived params (const)
    xi_th = sym('xt',[n,6]);
    xi_del = sym('xd',[n,6,6]);
    e_xi = sym('ex',[4,4,2*n]);
    K = sym('k',[6,6,n]);
    m = sym('m',[6,6,2*n]);
    m(:,:,:) = zeros(6,6,2*n);
    for i = 1:n
        %xi_theta
        xi_th(i,:) = [-cross(w_i(i,:),q_i(i,:)),w_i(i,:)];
        %xi_del
        xi_del(i,1:3,1:3) = eye(3);
        xi_del(i,4:6,4:6) = eye(3);
        xi_del(i,4:6,1:3) = zeros(3);
        xi_del(i,1:3,4) = -cross([1;0;0],q_i(i+1,:));
        xi_del(i,1:3,5) = -cross([0;1;0],q_i(i+1,:));
        xi_del(i,1:3,6) = -cross([0;0;1],q_i(i+1,:));
        %e_xi (run time, but done here for efficiency)
        e_xi(:,:, 2*i-1) = twistexp(xi_th(i,:).',q_th(i));
        e_xi(:,:,2*i) = eye(4);
        for k = twist_order
            e_xi(:,:,2*i) = e_xi(:,:,2*i)*twistexp(xi_del(i,:,k).',q_del(i,k));
        end
        %K
        K(:,:,i) = [E*A/l(i) 0 0 0 0 0;
                    0 3*E*IA_y/l(i)^3 0 0 0 (2*E*IA_y)/l(i)^2;
                    0 0 (3*E*IA_z)/l(i)^3 0 (2*E*IA_z)/l(i)^2 0;
                    0 0 0 (G*IA_r)/l(i) 0 0;
                    0 0 (2*E*IA_y)/l(i)^2 0 (E*IA_y)/l(i) 0;
                    0 (2*E*IA_z)/l(i)^2 0 0 0 (E*IA_z)/l(i)];
       
       %K(:,:,i) = [E*A/l(i) 0 0 0 0 0;
       %             0 3*E*IA_y/l(i)^3 0 0 0 0;
       %             0 0 (3*E*IA_z)/l(i)^3 0 0 0;
       %             0 0 0 (G*IA_r)/l(i) 0 0;
       %             0 0 0 0 (E*IA_y)/l(i) 0;
       %             0 0 0 0 0 (E*IA_z)/l(i)];
                
        m(1:3,1:3,2*i-1) = m_i(2*i-1)*eye(3);
        m(4:6,4:6,2*i-1) = I_i(:,:,2*i-1);
        m(1:3,1:3,2*i) = m_i(2*i)*eye(3);
        m(4:6,4:6,2*i) = I_i(:,:,2*i);
    end
%run-time params
    %sub and simplify first for computation time
    
    %need to update funs if you change these
    l_sub = [1.5,1.5].'; %meter
    r_sub = [.75, 0, 0;
            .75,0,0].';
    m_i_sub = [1,.1,1,2].'; %kgs
    I_sub = repmat(eye(3),[1,1,2*n]);
    I_sub(:,:,ind_del) = I_sub(:,:,ind_del)*.1;
    tan_del_sub = .05;
    m = subs(m,m_i,m_i_sub);
    m = subs(m,I_i,I_sub);
    
    diameter = .05;
    thickness = .005;
    K_params = [70e9,27e9, pi/4*((diameter)^2 - (diameter-2*thickness)^2),pi/64*diameter^4,pi/64*diameter^4, pi/32*diameter^4];
    K = subs(K,l,l_sub);
    K = subs(K,[ E G A IA_y IA_z IA_r],K_params);
    
    xi_th = subs(xi_th,l,l_sub);
    xi_th = double(subs(xi_th,r,r_sub));
    xi_del = subs(xi_del,l,l_sub);
    xi_del = double(subs(xi_del,r,r_sub));
    e_xi = subs(e_xi,l,l_sub);
    e_xi = simplify(subs(e_xi,r,r_sub));
    g_0 = subs(g_0,l,l_sub);
    g_0 = double(subs(g_0,r,r_sub));
        %g_si
    
    "Calculating FK"
    [g_si,g_si_th,g_si_del] = calculateFK(e_xi,g_0,n,ind_th,ind_del);
        %J_si
    "Calculating Jacobian"
    J_si = calculateJacobian(e_xi,xi_th,xi_del,g_0,n,true);
    %}
    %%calculating M, C, N%%
    Q_sym = Qfromq(q_th,q_del,dq_th,dq_del,n,true);
    
     %N
     
    [N, N_Y, N_dY] = calculatePotentialMat(g_si,J_si,q_th,q_del,m_i_sub,K,n);
    %M,C
    M = calculateMassMat(J_si,m,n);
    
    %assumption: bending does not effect coriolis ie, dMdq_del = 0
    [C, C_Y] = calculateCoriolisMat(M,q_th,dq_th,Q_sym,n,tan_del_sub);
   
    
    M_fun = createMassFuns(M,n,q_th,q_del,true);
    C_fun = createCoriolisFuns(C,n,q_th,q_del,dq_th,dq_del,false);
    %}
    N_fun = createPotentialFuns(N,n,q_th,q_del,true);    
    
    %[M_fun,C_fun,N_fun,T_fun] = getFunsFromFile(n);
    
    
    T = sym('T',[2*n,6]);
    target_pose = [pi/4,pi/4];
    F_ext = [0,0];
    T_th = sym('Tt',[n,1]);
    T(:,2:6) = zeros(2*n,5);
    T(ind_del,:) = zeros(n,6);
    T_fun = matlabFunction(T,'Vars',{T_th},'File','DynamicMats\T_fun')
    
    %{
    Y = C_Y;
    Y(1:7*n,1:7*n) = Y(1:7*n,1:7*n) + N_dY;%flattenDynamicMat(permute(C,[3,4,1,2]),n) + flattenDynamicMat(permute(Ydy_c,[3,4,1,2]),n) + N_dY;
    Y(1:7*n,7*n+1:end) = Y(1:7*n,7*n+1:end) + N_Y;%flattenDynamicMat(permute(C,[3,4,1,2]),n) + flattenDynamicMat(permute(Ydy_c,[3,4,1,2]),n) + N_dY;
    Y(7*n+1:end,1:7*n) = eye(7*n);
    Y(7*n+1:end,7*n+1:end) = zeros(7*n);
    
    "computing Y"
    Y_fun = matlabFunction(Y,'Vars',{Q_sym});
    %}

    %Y_ = jacobian([flattenDynamicMat(permute(C,[3,4,1,2]),n,true)*Q_sym(1:7*n)+flattenStateMat(N,n,true);Q_sym(1:7*n)],Q_sym);
    %Y_fun = matlabFunction(Y_,'Vars',{Q_sym})
    %Y_fun_mat = arrayfun(@(S)matlabFunction(S,'Vars',{Q_sym}),Y,'UniformOutput',false);
        
    th_i = ones(n,1).* (pi/6);
    del_i = [0,.05,0,0,0,0;0,0,0,0,0,0];
    dth_i = zeros(n,1);
    ddel_i = zeros(n,6);
    Q_init = Qfromq(th_i,del_i,dth_i,ddel_i,n)
    %[ E G A IA_y IA_z IA_r]
    %Using Aluminum 6061, 1m x .05m, t=.0033m in hollow cylinder 
        %calculates states in freefall simulation
     %evaluates g_si (for drawing purposes only)
    g_si_fun = matlabFunction(subs(subs(g_si,l,l_sub),r,r_sub), 'vars', {q_th,q_del});
    %joints and CoM of the robot
    J_si_fun = matlabFunction(subs(subs(J_si,l,l_sub),r,r_sub), 'vars', {q_th,q_del});
    tspan = 0:.001:10;
    figure(2)
    clf(2)
    figure(3)
    clf(3)
    figure(4)
    clf(4)
    fcelleval = @(F,args)feval(F{1},args{:});
    %displayState(Q_init,g_si_fun,J_si_fun,n,r_sub,ind_th,ind_del,'Jacobian',true,"")
    
   
%{
     ranges = [linspace(0,2*pi,40);
              linspace(-.05, .05,40);
              linspace(-.05, .05,40);
              linspace(-.05, .05,40);
              linspace(0,2*pi,40);
              linspace(0,2*pi,40);
              linspace(0,2*pi,40)]
          for i = 1:2*n
        for j = 1:7
            q_draw = zeros(14*n,1);
            for r = ranges(j,:)
                [i,j-1,r]
                q_draw(7*i-7+j) = r;
                displayState(0,q_draw,g_si_fun,J_si_fun,n,r_sub,ind_th,ind_del,'velocity',false,'')
            end
        end
    end
  %}  
    [t_freefall,y_freefall] = ode15s(@(t,y)EoM(t,y,C_fun,N_fun,@(qt,qd,dqt,dqd) zeros(2*n,6),n,ind_th,ind_del,double(K)),...
                                tspan,...
                                Q_init,...
                                odeset('Mass',@(t,y)Mass_eval(t,y,M_fun,n,ind_th,ind_del),...
                                       'OutputFcn',@(t,Q,msg)displayState(t,Q,g_si_fun,J_si_fun,n,r_sub,ind_th,ind_del,'velocity',false,msg)));%,...
                                       %'Jacobian',@(t,Q)Y_fun(Q)));
                                       %'Jacobian',@(t,Q)arrayfun(@(F)fcelleval(F,{Q}),Y_fun_mat)));


    p = zeros(4,2*n,size(y_freefall,1));
    for sim = 1:1
        if sim == 1
            y = y_freefall;
            t = t_freefall;
            "running freefall simulation"
            %pause()
        else
            %y = y_controlled;
            %t = t_controlled;
            %"running controller simulation"
            %pause()
        end
   
        for t_ = 1:size(y,1)
            figure(1)
            hold on
            grid on
            time = t(t_)
            displayState(y(t_),g_si,J_si,n,r_sub,ind_th,ind_del,'none',false)
        end
    end
   
  th_s = [pi/6,pi/6].';
  del_s = [.1,0,0,0,0,0;
            0,0,.1,0,0,0];
        
  function status = displayState(t, Q, g_si_fun,J_si_fun,n,r,ind_th,ind_del,vecshow, midpause, msg)
    if nargin == 8
        vecshow='none';
        midpause = false;
    elseif nargin==9
        midpause=false;
    elseif nargin==10
        msg='';
    end
    if strcmp(msg,'init')
        status = 0;
        return
    end
    if strcmp(msg,'init')
        return
        vecshow='Jacobian';
    elseif strcmp(msg,'done')
    end
    [q_s,dq_s,th_s, del_s, dth_s,ddel_s] = Qtoq(Q(:,end),n,ind_th,ind_del);
    plt_del1 = figure(3);
    t(1)
    hold on
    scatter(repmat(t(1),1,6),del_s(1,:),10,[255 0 0;0 255 0;0 0 255;255 255 0; 255 0 255; 0 255 255])
    hold off
    plt_th = figure(4);
    hold on
    scatter(repmat(t(1),1,2),th_s.',10,[255 0 255; 0 255 255])
    hold off
    robot = figure(2);
    %l_ = [0,5,5].';
    g_si_q0 = g_si_fun(th_s,del_s);
    J_si_q0 = J_si_fun(th_s,del_s);
    g_base = zeros(4,2*n);
    g_end = zeros(4,2*n);
    for i = 1:n
        g_base(:,2*i-1) = [-r(:,i);1];
        g_base(:,2*i) = [0,0,0,1].';
        g_end(:,2*i-1) = [r(:,i);1];
        g_end(:,2*i) = [0,0,0,1].';
    end    
    if strcmp(vecshow,'Jacobian')
        vecs_th = zeros(4,2*n,n); %vec_len, vecs/set, vec sets
        vecs_del = zeros(4,6*n,n);
        vecs_th(1:3,:,:) = permute(J_si_q0(1:3,ind_th,1,:),[1,4,2,3]);
        vecs_del(1:3,:,:) = reshape(permute(J_si_q0(1:3,ind_del,1:3,:),[1,3,4,2]),3,6*n,n);
    elseif strcmp(vecshow,'velocity')

        vels = squeeze(sum(dmatmul(dq_s,J_si_q0,true,[3,2],[3,1]),3));
    end
    p = zeros(4,3*2*n);
    for i = 1:2*n
        p(:,3*i-2) = g_si_q0(:,:,i)*g_base(:,i);        
        p(:,3*i-1) = g_si_q0(:,:,i)*[0,0,0,1].';
        p(:,3*i) = g_si_q0(:,:,i)*g_end(:,i);
    end
    p = p(1:3,:);
    p_mid = p(:,mod(1:3*n*n,3)==2);
    clf(robot)
    plot3(p(1,:),p(3,:),p(2,:),'-o','color','blue')
    hold on
    if strcmp(vecshow,'Jacobian')
        vect_fig = gobjects(2*n)
        for i = 1:n
            for j = 1:2*n
                vecs_th(:,j,i) = g_si_q0(:,:,j) * vecs_th(:,j,i);
            end
            vect_fig(2*i-1) = quiver3(p_mid(1,:),p_mid(2,:),p_mid(3,:),vecs_th(1,:,i),vecs_th(2,:,i),vecs_th(3,:,i));
            for j = 1:2*n
                vecs_del(:,3*j-2,i) = g_si_q0(:,:,j)*vecs_del(:,3*j-2,i);
                vecs_del(:,3*j-1,i) = g_si_q0(:,:,j)*vecs_del(:,3*j-1,i);
                vecs_del(:,3*j,i) = g_si_q0(:,:,j)*vecs_del(:,3*j,i);
            end
            p3 = repelem(p_mid,1,3);
            vect_fig(2*i) = quiver3(p3(1,:),p3(3,:),p3(2,:),vecs_del(1,:,i),vecs_del(3,:,i),vecs_del(2,:,i));
        end
    elseif strcmp(vecshow,'velocity')
        vel_vecs = zeros(4,4,2*n);
        for i = 1:2*n
            V = [auto_skew(vels(4:6,i)), vels(1:3,i);0,0,0,0];
            vel_vecs(:,:,i) = g_si_q0(:,:,i) * V;
        end
        scale = max(max(squeeze(abs(vel_vecs(1:3,4,:)))));
        if scale ~= 0
            scale = 1/scale;
        end
        vecs = squeeze(vel_vecs(1:3,4,:))*scale;
        vect_fig = quiver3(p_mid(1,:),p_mid(3,:),p_mid(2,:),vecs(1,:),vecs(3,:),vecs(2,:),'color','red');
        
    end
    zlim([-3 4])
    xlim([-4,4])
    ylim([-1,1])
    view(-22,51)
    box on
    view(-7,5)
    hold off
    if vecshow
        dcm_obj = datacursormode(robot);
        map = repelem(1:2*n,3);
        set(dcm_obj,'UpdateFcn',{@poll_joint,robot,vect_fig,map})
    end
    if midpause
        pause()
    end
    drawnow
    status = 0;
  end
  