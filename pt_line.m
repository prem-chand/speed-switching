function [obj]=pt_line(pt0,pt1,color,thickness,parent)
obj=line([pt0(1) pt1(1)],[pt0(2) pt1(2)],'Parent',parent);
set(obj,'Color',color);
set(obj,'LineWidth', thickness);