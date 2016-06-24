%function sliding_T
% MFPT of a sliding protein. means ksize=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% there are no traps in the path: add=1 or nt=1 and xp=xA;

%     ******          T = N^2 + 2(N-xp)(tau-1)      ******

xZ=0; xA=25; N=25;
ksize=1;
add=1;

for nt=2:5
    file=sprintf('TeffectR_xZ%d_xA%d_N%d_nt%d_add%d_ksize%d.txt',xZ,xA,N,nt,add,ksize);
    if(exist(file) > 1)
        a=importdata(file);
        aval(nt)=mean(a(:,2));
    else
        disp(file);disp('does not exist');
    end
end

%% effect of single trap with different depths.

%     ******          T = N^2 + 2(N-xp)(tau-1)      ******

clear;
clc;

xZ=0; xA=25; N=25;
ksize=1;

nt=1;
i=1;
for add=1:10
    file=sprintf('TeffectR_xZ%d_xA%d_N%d_nt%d_add%d.txt',xZ,xA,N,nt,add);
    if(exist(file) > 1)
        a=importdata(file);
        x=a(2:end,1);
        y=a(2:end,2);
        [g,h]=fit(x,y,'Poly1');
        aval(i)=g.p1;
        bval(i)=g.p2;
        i=i+1;
    else
        disp(file);disp('does not exist');
    end
end


%% traps (multi) with different energy depths
%for each nt the fit will vary since the effect of trap positions:
% Link to GapNDistanceOfTrappNT(nt).m