%%%%%%%%%%%%%%%%%%%%%%%%%%
% finding the percentage of proteins that has highest trap distribution
% around target site.
% inputs : trap positions file
%%%%%%%%%%%%%%%%%%%%%%%%%%

%function maxTrap_prot_percentage
clear;
clc;

sp='mouse';
mats=127;

around=1000;
s_tot=around/100;

path=sprintf('/Volumes/Jeni_Seagate/TRAP-11/%s',sp);
addpath(path);
%mat = importdata(sprintf('%s_matno.txt',sp));
%mats=size(mat,1);
j=1;
X=(0-around-50):100:(around-50);

for i=1:mats
    %mm=mat(i,1);
    mm=i;
    file=sprintf('Trap_%s_%d.txt',sp,mm);
    if(exist(file) > 0)
        a=importdata(file);
        a(a==0)=[];

        if(isempty(a) == 0)
            [n,out]=hist(a,X);
            [~,m]=max(n);
            b(j,1)=out(m);
            j=j+1;
        end
    end
end

%************************** Percentage for each gap
c=zeros(s_tot,1);
for i=1:s_tot
    min_r=((i-1)*100)+1;
    max_r=(i*100);
    max_l=0-min_r;
    min_l=0-max_r;
    
    for j=1:size(b,1)
        if((b(j,1) >=min_r && b(j,1) <= max_r) || (b(j,1) >=min_l && b(j,1) <= max_l))
            c(i,1)=c(i,1)+1;
        end
    end
end
s=sum(c);
c=c/s;
fm=fopen(sprintf('maxTrap_pro_dist_%s.txt',sp),'wt');
for i=1:s_tot
    fprintf(fm,'%s,%0.3f\n',sprintf('around %d',(i*100)),c(i,1));
end
fclose(fm);
