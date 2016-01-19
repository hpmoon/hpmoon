function [v, v1]=hypervolume(P1,R)
P = sortrows(P1,1);
c = size (P,1);
PT (1,:) = R(1,:);
PT(2,:) = [P(1,1),0];
PT = vertcat (PT,P);
c1 = size (PT,1);
PT (c1+1,1)= R(1,1);
PT (c1+1,2)= P(c,2);
PT (c1+2,:) = R(1,:); 
base = 0;
ar = polyarea(PT(:, 1), PT(:, 2));
%figure;
%plot(PT(:, 1), PT(:, 2));
for i=1 : c
    if i == 1
    d1= abs((R(1,1)-P(i,1)));
    d2= abs((R(1,2)-P(i,2)));
    base = base + (d1*d2);
    else
    d1= abs((R(1,1)-P(i,1)));
    d2= ((P(i-1,2)-P(i,2)));
    base = base + (d1*d2);
    end;
end
v=ar;
v1 = base;