GRF_h=GRF(:,1);
GRF_v=GRF(:,2);
hold on;
for i=1:1:length(GRF_h)
    mag=norm(GRF(i,:));
    plot([0 GRF_h(i)/mag],[0 GRF_v(i)/mag]);
    plot(cm(i,1),cm(i,2),'x');
end
hold off;