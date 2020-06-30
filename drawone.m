function drawone(p,parent,shift_vec)

tem = get(parent,'Children');

delete(tem);

pstk    = p(1:2,1);
ph      = p(3:4,1);
pswk    = p(5:6,1);
pswf    = p(7:8,1);
phead   = p(9:10,1);

%[pF_L,pK_L,pF_R,pK_R,pHead,pThigh_COM_L,pShin_COM_L,pThigh_COM_R,pShin_COM_R,pTorso_COM,qThigh_L,qShin_L,qThigh_R,qShin_R]=points(q,pH,[]);

%[pF,pK,pHead,pThigh_COM,pShin_COM,pTorso_COM,qThigh,qShin]=points(q,pH,[]);

% BiMASC configuration - Stance only!
% pF    = [0;0]; % how can we include the flight??
% pH    = fcn_stance_pHip(q);
% pK    = fcn_stance_pKnee(q);
% pHead = fcn_stance_pHead(q);
% 

gray = [ 0.500 0.500 0.500 ];


shift_vec1 = shift_vec';%[shift_vec;0];


pstk = pstk + shift_vec1;

ph   = ph + shift_vec1;

pswk = pswk + shift_vec1;

pswf = pswf + shift_vec1;

phead= phead + shift_vec1;


pt_link(shift_vec1,pstk,'g','g',parent);

pt_link(pstk,ph,'g','g',parent);

pt_link(ph,pswk,'g','g',parent);

pt_link(pswk,pswf,'g','g',parent);

pt_link(ph,phead,'b','b',parent);


radius = 0.04;

pt_circle(pstk,radius,'k',parent);
pt_circle(pswk,radius,'k',parent);
pt_circle(ph,radius,'k',parent);
pt_circle(phead,0.06,'k',parent);

%pt_line([-10 0],[10 0],'k',0.5,parent);
o1=[0.4 0];
pt_line([-10 0],o1,'k',0.5,parent);
for i=1:5
    o2=o1+[0 0.1];
    pt_line(o1,o2,'k',0.5,parent);
    o3=o2+[0.22 0];
    pt_line(o2,o3,'k',0.5,parent);
    o1=o3;
end
pt_line(o1,o1+[10 0],'k',0.5,parent);



%pt_circle(pF,radius,'c',parent);