%*************************************
%function TrapVsBindsize
% compare the number of traps of a matfile with its bind site size
% infile : Bindsize_species.txt => matno,size of bindsite, no of traps
% *************************************
clear;
clc;
clf;

sp=importdata('species.txt');

for s=1:size(sp,1)
    outt=[];
    species=sp{s,1};
    path=sprintf('/Volumes/Jeni_Seagate/TRAP-11/%s',species);
    addpath(path);

    infile_up=sprintf('bindsizeREV_traps_%s_up.txt',species);
    infile_down=sprintf('bindsizeREV_traps_%s_down.txt',species);

    outfile=sprintf('hist_bsize_trap_%s.txt',species);

    a=importdata(infile_up);
    a1=importdata(infile_down);
    a=[a;a1];
    A=a(:,3);
    A(A==-1)=0;
    a(:,3)=A;
    
    rows = size(a,1);

    str=min(a(:,2));
    ed=max(a(:,2));
    out=zeros(ed-str+1,1);
        
    for j=1:rows
        for i=str:ed
            if(a(j,2) == i)
                out(i-str+1,1)=out(i-str+1,1)+a(j,3);
                break;
            end
        end
    end
    
    aa=[];
    summ = sum(out);
    aa = out/summ;
    out=aa;
    outt(:,1)=str:ed;
    outt(:,2)=out;
    
    outfile=sprintf('TrapsDist_%s.txt',species);
    dlmwrite(outfile,outt);
end

% fig=figure;
% axes1=axes('Parent',fig);
% box(axes1,'on');
% hold(axes1,'all');
% bar3(axes1,[1:maxx],out(1:maxx,:));
% 
% title('3000 around TSS','FontSize',14);
% 
% legend(sp);