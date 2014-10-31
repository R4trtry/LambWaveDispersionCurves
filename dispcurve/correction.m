
% fix the curve by sorting and matching the order
for j = Lf-1:-1:1    
    for i = acount:-1:1
        [~,I]=min(abs(dispca(j,i)-dispca(j+1,:)));
        [~,I1] = min(abs(dispca(j+1,I) - dispca(j,:)));
        if(~isnan(I) && I~=i && I1==i), 
            dispca(j,I) = dispca(j,i); dispca(j,i) = nan; 
        end
    end 
    for i = scount:-1:1
        [~,I]=min(abs(dispcs(j,i)-dispcs(j+1,:)));
        [~,I1] = min(abs(dispcs(j+1,I) - dispcs(j,:)));
        if(~isnan(I) && I~=i && I1==i), 
            dispcs(j,I) = dispcs(j,i); dispcs(j,i) = nan; 
        end
    end
end
