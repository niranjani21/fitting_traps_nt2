function NT2_plotting
% needs saved coefficient file

%% First coefficient 
a=importdata('coeff1_nt2.txt');
xp1=2; 
k=1; j=2;

for i=xp1+1:20
    val(i-xp1,k)=(a(j,2)*xp1)+(a(j,3)*i)+a(j,4);
end

x=3:20;
y=val;
g=polyfit(x',y,1);
yy=polyval(g,x');
plot(x',yy,'LineWidth',2.0); hold on;
plot(x',y,'*','MarkerSize',10.0);
xlabel('Trap position ($$ x_{p^3} $$)','Interpreter','Latex',...
    'FontName','Times New Roman','FontSize',12,'FontWeight','b');
ylabel('MFPT','Interpreter','Latex','FontName','Times New Roman',...
    'FontSize',12,'FontWeight','b');
h=text('Interpreter','Latex','String','$$ T = ax_{p^1}+bx_{p^2}+c $$','Position',[11 715],'FontSize',13,'FontWeight','b');
print -n painters -dpng -r600 nt2_coeff1.png

%*******************************************
%% second coefficient
c=importdata('coeff1_nt2.txt');

subplot(2,2,1); 
x=c(:,1);
y=c(:,2);
g=polyfit(x,y,1);
yy=polyval(g,x);
plot(x,yy,'r','LineWidth',2.0); hold on;
plot(x,y,'*','MarkerSize',10.0);
xlabel('a','Interpreter','Latex','FontName','Times New Roman','FontSize',10,'FontWeight','b');

subplot(2,2,2);
y=c(:,3);
g=polyfit(x,y,1);
yy=polyval(g,x);
plot(x,yy,'r','LineWidth',2.0); hold on;
plot(x,y,'*','MarkerSize',10.0);
xlabel('b','Interpreter','Latex','FontName','Times New Roman','FontSize',10,'FontWeight','b');

subplot(2,2,[3 4]);
y=c(:,4);
g=polyfit(x,y,1);
yy=polyval(g,x);
plot(x,yy,'r','LineWidth',2.0); hold on;
plot(x,y,'*','MarkerSize',10.0);
xlabel('c','Interpreter','Latex','FontName','Times New Roman','FontSize',10,'FontWeight','b');

ylabel('$$ \tau_t -1 $$','Interpreter','Latex','FontName','Times New Roman','FontSize',10,'FontWeight','b');
title('Inset I','Interpreter','Latex','FontName','Times New Roman','FontSize',10,'FontWeight','b');
print -n painters -dpng -r600 nt2_coeff2.png

%******************************************
%%
a=importdata('coeff1_nt2.txt');
xp1=2; 
k=1;
for j=3:2:10
    txt='\tau_t';
    leg{k,1}=sprintf('%s : %d',txt,j);
    for i=xp1+1:20
        val(i-xp1,k)=(a(j,2)*xp1)+(a(j,3)*i)+a(j,4);
    end
    k=k+1;
end
plot(val,'LineWidth',2.0); legend(leg);
xlabel('Trap position ($$ x_{p^3} $$)','Interpreter','Latex',...
    'FontName','Times New Roman','FontSize',12,'FontWeight','b');
ylabel('MFPT','Interpreter','Latex','FontName','Times New Roman',...
    'FontSize',12,'FontWeight','b');
h=text('Interpreter','Latex','String','$$ T = ax_{p^1}+bx_{p^2}+cx_{p^3}+d $$','Position',[14 756],'FontSize',13,'FontWeight','b');
print -n painters -dpng -r600 nt2_coeff3.png
