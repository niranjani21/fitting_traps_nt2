function GapNDistanceOfTrappNT2(xZ,nt)
% gap between two traps and avg distance of traps from absorbing boundary
% inputs : number of traps, position of absorbing boundary.
%%%%%%%%%%%%%%%%

combfile=sprintf('comb_%d.txt',nt);
combin=importdata(combfile);

n=nt-1;
j=1;
for i=25:-1:n
    nu=factorial(i);
    de=factorial(n)*factorial(i-n);
    val=nu/de;
    if(i==25)
        com(j)=val;
        j=j+1;
    else
        com(j)=com(j-1)+val;
        j=j+1;
    end
end
rows=size(com,2);

for ksize=1:8
    coval2=[];
    for add=1:10
        file=sprintf('TeffectR_xZ%d_xA25_N25_nt%d_add%d_ksize%d.txt',xZ,nt,add,ksize);
        a=importdata(file);
        coval=[];
        if(add==1)
            coval2(add,3)=mean(a(:,2));
        else
            
            for i=2:20              % i starts at 2(1) cause, 1(0) is the initial position of the protein.
                x=[]; y=[];
                if(i==1)
                    st=1;
                else
                    st=com(i-1)+1;
                end
                ed=com(i);
                x=combin(st:ed,nt);
                y=a(st:ed,2);
                [g,~]=fit(x,y,'Poly1');
             
                coval(i-1,1)=g.p1;
                coval(i-1,2)=g.p2;

            end
            coval2(add,1)=mean(coval(:,1));
            X=2:20;
            [g1,~]=fit(X',coval(:,2),'Poly1');
            coval2(add,2)=g1.p1;
            coval2(add,3)=g1.p2;
        end
    end
    A=[0:9];
    [k1,~]=fit(A',coval2(:,1),'a*x');
    [k2,~]=fit(A',coval2(:,2),'a*x');
    [k3,~]=fit(A',coval2(:,3),'Poly1');
    
    y1(ksize,1)=ksize;
    y1(ksize,2)=k1.a;
    y1(ksize,3)=k2.a;
    y1(ksize,4)=k3.p1;
    y1(ksize,5)=k3.p2;
end
K=[1:8];
[n1,hn1]=fit(K',vals(:,1),'(6*a)/((x+1)*(2*x+1))');
[n2,hn2]=fit(K',vals(:,2),'(6*a)/((x+1)*(2*x+1))');
[n3,hn3]=fit(K',vals(:,3),'(6*a)/((x+1)*(2*x+1))');
[n4,hn4]=fit(K',vals(:,4),'(6*625)/((x+1)*(2*x+1))+(b/x)+((3*x*(x+1))/(2*(2*x+1)))');

out(1,1)=n1.a;
out(1,2)=n2.a;
out(1,3)=n3.a;
out(1,4)=n4.a;
out(1,5)=n4.b;
 
% outfile=sprintf('coeff_nt%d.txt',nt);
% dlmwrite(outfile,y1);


%% plotting

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

%%
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

