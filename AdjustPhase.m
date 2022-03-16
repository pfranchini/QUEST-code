phaseadjust =1;
str = clipboard('paste');
%str = sprintf ('4 85 7\n4 7 8\n');

data = textscan (str, '%f\t%f\t%f\n');
freq = data{1}';
X = 1e3*data{2}';
Y = (1e3*data{3}');

cut=2;
cutend=12;
X=X(cut:end-cutend);
Y=Y(cut:end-cutend);
freq = freq (cut:end-cutend);

figure(1)
clf
hold on
p1=plot(freq,X,'ob');
p2=plot(freq,Y,'or');


% Define boundaries for fit parameters    
    lb = [-1e3,-1e3,0,0,0000,-phaseadjust*pi,-1e3,-1e3];
    ub = [1e3,1e3,inf,1e7,11e7,phaseadjust*pi,1e3,1e3];
% Fit options    
    options = optimoptions('lsqcurvefit','TolFun',1e-24,'TolX',1e-24,'MaxFunEvals',1e4,...
        'MaxIter',1e2);
% f=1e6;
% d=1e4;
    aa = [];
    msee = [];

% Try fits with different initial phases b
    for k = 1:8%8;

        b = -pi+k*pi/4;
   
        X1 = X*cos(b) + Y*sin(b);
        Y1 = -X*sin(b) + Y*cos(b);
     % Initial Amplitude
    [A,ind] = max(X1);
    A = (A - min(X1)); 

% Initial Width
    [B,ind1] = max(Y1);
    [B,ind2] = min(Y1);
    f = freq(ind);
    d = abs(freq(ind1)-freq(ind2))/8;   
       
  
    a0=[mean(X1),0,A,d,f,b,mean(Y1),0];
   % a0 = [0,0,10, 50, 16500, b, 0,0];
    [a, mse] = lsqcurvefit(@lorv,a0,freq,[X;Y],lb,ub,options);
    mse;
    aa=[aa;a];
    msee=[msee;mse];
    end;

% Which fit gives smallest error?
    [mse, ind] = min(msee);
   
    a = aa (ind,:);
% Output curve to clipboard
% T = table(freq',X',Y');
% TT = sprintf('%f\t%f\t%f\r\n',[freq;X;Y]);
% clipboard ('copy',TT);

% Phase in degrees from -180 to +180

f0 = a(5)
width = a(4)
height = a(3) *1e3
hw =  width * height 
phase = mod((a(6))/pi*180+180,360)-180
% plot best fit

F_fit = lorv(a,freq);

p3=plot(freq, F_fit(1,:),'b','LineWidth',2);
p4=plot(freq,  F_fit(2,:),'r','LineWidth',2);

legend([p1,p2],'In-phase voltage','Out-of-phase voltage')
    

    
    
    
    
    
    
