function [newlist] = kat_stats_commonality_genlist(ivlist,value)

numlist = length(ivlist);
newlist = ivlist;
newlist = 0;
for i = 1:numlist
    newlist(i) = abs(ivlist(i)) + abs(value);
    if (((ivlist(i) < 0) && (value >= 0)) || ((ivlist(i) >=  0) && (value < 0))) 
        newlist(i) = newlist(i) * -1;
    end
end
