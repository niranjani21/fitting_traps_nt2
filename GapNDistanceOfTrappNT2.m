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
[n4,hn4]=fit(K',vals(:,4),'(6*625)/((x+1)*(2*x+1))+(b/x)+((3*x*(x+1))/(2*(2*x+1)))');

out(1,1)=n4.a;
out(1,2)=n4.b;
 
% outfile=sprintf('coeff_nt%d.txt',nt);
% dlmwrite(outfile,y1);
