%----------------------------------------------------------------------------%
% AS5430 COURSE-PROJECT                                                      %
% LINEAR INVISCID INSTABILITY STUDY OF TANH VELOCITY PROFILE                 %
% GNU OCTAVE / MATLAB POST-PROCESSING FILE                                   %
%----------------------------------------------------------------------------%

%--------plot 1: temporal growth rate vs wavenumber--------------------------%
% data:
file1  = '/home/vikas/Desktop/AS_5430_project_ME15D080/output/disp_data';
disp_data = load(file1,'-ascii');
k  = [0; disp_data(:,1);1]; K = linspace(0,1); % wavenumber
g  = [0; disp_data(:,2);0]; G = spline(k,g,K); % temporal growth

% figure:
f1 = figure(1);
W = 4; H = 4;
set(f1,'PaperUnits','inches');set(f1,'PaperOrientation','portrait');
set(f1,'PaperSize',[H,W])    ;set(f1,'PaperPosition',[0,0,W,H]);
plot(K,G,'k-','LineWidth',1.5);axis([0,1,0,0.1]);
xlabel('\kappa');ylabel('\kappa c_{i}');
print(f1,'-deps','-color','plot 1: disp_data.eps');


%--------plot 2 & 3: phi_eigen_fcns vs y ------------------------------------%
% data
for count = 2:2:8
  index=count/2;
  temp_file  = sprintf('/home/vikas/Desktop/AS_5430_project_ME15D080/output/%d/eig_fcn', count);
  PHI_EIG_DATA {index}=load(temp_file,'-ascii');
end

% figures:
%--------------------plot -2 real_phi ---------------------------------------%
f2 = figure(2);
W = 5; H = 3;
set(f2,'PaperUnits','inches');set(f2,'PaperOrientation','portrait');
set(f2,'PaperSize',[H,W])    ;set(f2,'PaperPosition',[0,0,W,H]);

for i= 1:4
 eig_data=PHI_EIG_DATA{i};
 plot(eig_data(:,1),eig_data(:,2),'k-','LineWidth',1.5);axis([0,7,0,1.2]);
 xlabel('y');ylabel('\phi_{r}');
 print(f2,'-deps','-color','plot 2: phi_r.eps');
 hold on;
end
 hold off;

%--------------------plot -3 imag_phi----------------------------------------%
f3 = figure(3);
W = 5; H = 3;
set(f3,'PaperUnits','inches');set(f3,'PaperOrientation','portrait');
set(f3,'PaperSize',[H,W])    ;set(f3,'PaperPosition',[0,0,W,H]);

for i= 1:4
 eig_data=PHI_EIG_DATA{i};
 plot(eig_data(:,1),eig_data(:,3),'k-','LineWidth',1.5);axis([0,7,0,1.2]);
 xlabel('y');ylabel('\phi_{i}');
 print(f3,'-deps','-color','plot 3: phi_i.eps');
 hold on;

end

 hold off;
%-----------------plot -4 & 5 omg_eigen_fcns vs y----------------------------%
% data
iota=sqrt(-1);
for i = 1:4
  eig_data = PHI_EIG_DATA{i};
  lambda   = (tanh(eig_data(:,1)).*(sech(eig_data(:,1))).^2)./(0.5*tanh(eig_data(:,1))-iota*disp_data(2*i,3));
  OMG_DATA{i} = lambda.*(eig_data(:,2)+iota*eig_data(:,3));
end
%------------------ plot-4 real_omg -----------------------------------------%
f4 = figure(4);
W = 5; H = 3;
set(f4,'PaperUnits','inches');set(f4,'PaperOrientation','portrait');
set(f4,'PaperSize',[H,W])    ;set(f4,'PaperPosition',[0,0,W,H]);

for i= 1:4
 omg_data= OMG_DATA{i};
 plot(eig_data(:,1),real(omg_data),'k-','LineWidth',1.5);axis([0,1.8,0,2]);
 xlabel('y');ylabel('\omega_{r}');
 print(f4,'-deps','-color','plot 4: omg_r.eps');
 hold on;
end
 hold off;
%------------------ plot-5 imag_omg ----------------------------------------%
f5 = figure(5);
W = 5; H = 3;
set(f5,'PaperUnits','inches');set(f5,'PaperOrientation','portrait');
set(f5,'PaperSize',[H,W])    ;set(f5,'PaperPosition',[0,0,W,H]);

for i= 1:4
 omg_data= OMG_DATA{i};
 plot(eig_data(:,1),imag(omg_data),'k-','LineWidth',1.5);axis([0,1.8,0,1.2]);
 xlabel('y');ylabel('\omega_{i}');
 print(f5,'-deps','-color','plot 5: omg_i.eps');
 hold on;
end
 hold off;
%--------------------screen display-----------------------------------------%
% DATA FOR MAX AMPLIFIED WAVE
 wave_number = disp_data(:,1);ww=linspace(disp_data(1,1),disp_data(end,1));
 temp_growth = disp_data(:,2);gg=spline(wave_number,temp_growth,ww);
 init_cond   = disp_data(:,3);ic=spline(wave_number,init_cond,ww);
 
 [val,ind] = max(gg);
 disp('%------------------------------------------------------------------%')
 disp('MOST AMPLIFIED WAVENUMBER = ');disp(ww(ind));
 disp('MAX GROWTH RATE = ')          ;disp(gg(ind));
 disp('CORRESPONDING Ci= ')          ;disp(gg(ind)/ww(ind));
 disp('%------------------------------------------------------------------%')
%----------plot 6 & 7 streamlines and vorticity contours--------------------%
% data
  max_file  = sprintf('/home/vikas/Desktop/AS_5430_project_ME15D080/k_max_data/k_max/eig_fcn');
  max_k     = load(max_file,'-ascii');
  max_k_grid= max_k;
  for i = 1:701
      max_k_grid(i,:)=max_k(702-i,:);
  end 

% make grid
  x =linspace(0,14,700);y =linspace(2.2,-5,701);
  [X,Y]=meshgrid(x,y);
  epsilon1=0.1;k0=0.4478;ci0=0.2120;
  lamb   = (tanh(max_k_grid(:,1)).*(sech(max_k_grid(:,1))).^2)./(0.5*tanh(max_k_grid(:,1))-iota*ci0);
  phi_r= zeros(701,700);phi_i= zeros(701,700);LAMB=zeros(701,700);
  for j=1:700
      phi_r(:,j)=max_k_grid(:,2);
      phi_i(:,j)=max_k_grid(:,3);
      LAMB (:,j)=lamb;
  end
  omg = LAMB.*(phi_r+iota*phi_i);
  str_fcn =   Y + 0.5*log(1+exp(-2*Y))+  epsilon1*(phi_r.*cos(k0*X)-phi_i.*sin(k0*X));
  vort_fcn= -0.5*sech(Y).*sech(Y)     +2*epsilon1*((real(omg).*cos(k0*X))-(imag(omg).*sin(k0*X)));

%------------------ plot-6 streamlines ----------------------------------------%
f6 = figure(6);
W = 5; H = 3;
set(f6,'PaperUnits','inches');set(f6,'PaperOrientation','portrait');
set(f6,'PaperSize',[H,W])    ;set(f6,'PaperPosition',[0,0,W,H]);
contour(X,Y,str_fcn,50,'LineWidth',1.5);axis([0,14,-5,0.5]);
xlabel('x');ylabel('y');
print(f6,'-deps','plot 6: phi_max.eps');

%------------------ plot-7 vorticity contours ----------------------------------------% 

f7 = figure(7);
W = 5; H = 3;
set(f7,'PaperUnits','inches');set(f7,'PaperOrientation','portrait');
set(f7,'PaperSize',[H,W])    ;set(f7,'PaperPosition',[0,0,W,H]);
contour(X,Y,vort_fcn,5,'LineWidth',1.5);axis([0,14,-2,2]);
xlabel('x');ylabel('y');
print(f7,'-deps','plot 7: vort_max.eps');

%--------------------------------------------------------------------------------------%
%  STREAMLINES AND VORTICITY CONTOURS FOR NEUTRAL CASE WITH KNOWN ANALYTIC SOLUTION    %
%--------------------------------------------------------------------------------------%

%------------------ plot-8 streamlines ----------------------------------------%
x =linspace(0,14);y =linspace(0.5,-5);
[X,Y]=meshgrid(x,y);eps_1=0.1;
str_fcn=Y+0.5*log(1+exp(-2*Y))+eps_1*(sech(Y).*cos(k0*X));

f8 = figure(8);
W = 5; H = 3;
set(f8,'PaperUnits','inches');set(f8,'PaperOrientation','portrait');
set(f8,'PaperSize',[H,W])    ;set(f8,'PaperPosition',[0,0,W,H]);
contour(X,Y,str_fcn,50,'LineWidth',1.5);axis([0,14,-5,0.5]);
xlabel('x');ylabel('y');
print(f8,'-deps','plot 6: phi_neutral.eps');

%------------------ plot-9 vorticity contours --------------------------------% 
x =linspace(0,14);y =linspace(2,-2);
[X,Y]=meshgrid(x,y);eps_2=0.2;
vort_fcn=-0.5*sech(Y).^2+ eps_2*(2*(sech(Y)).^3.*cos(X));
f9 = figure(9);
W = 3; H = 4;
set(f9,'PaperUnits','inches');set(f9,'PaperOrientation','portrait');
set(f9,'PaperSize',[H,W])    ;set(f9,'PaperPosition',[0,0,W,H]);
contour(X,Y,vort_fcn,5,'LineWidth',1.5);axis([0,6,-2,2]);
xlabel('x');ylabel('y');
print(f9,'-deps','plot 7: vort_neutral.eps');




