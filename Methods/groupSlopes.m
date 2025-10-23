function [labels, idx] = groupSlopes(slope_amp, options)

arguments
    slope_amp double
    options.methods (1,1) string = "kmeans"   % "kmeans"|"knn" or "pct"
    options.pct    (1,1) double {mustBeGreaterThan(options.pct,0),mustBeLessThan(options.pct,0.5)} = 0.2
end

A = slope_amp; flipBack = isvector(A) && iscolumn(A); if isvector(A), A=A(:)'; end
[M,T]=size(A); labels = nan(M,T);

for r=1:M
    s=A(r,:)'; v=~isnan(s);
    if any(v)
        switch lower(options.methods)
            case {"kmeans","knn"}
                if nnz(v)>=3
                    [c,mu]=kmeans(s(v),3,'Replicates',5,'MaxIter',100);
                    [~,ord]=sort(mu); map=-1:1; lbl=zeros(nnz(v),1);
                    for k=1:3, lbl(c==ord(k))=map(k); end
                    labels(r,v)=lbl;
                end
            otherwise  % "pct": tails = top/bottom pct; stable = pct closest to 0 from the remainder
                sv = s(v); n = numel(sv);
                k = max(1, round(options.pct*n));          % target count per group
                kTB = min(k, floor(n/2));                  % cap tails if n is small

                [~,ord] = sort(sv,'ascend');
                dec = ord(1:kTB);                          % bottom pct
                inc = ord(end-kTB+1:end);                  % top pct

                rem = true(n,1); rem([dec;inc]) = false;   % remove tails before picking stable
                remIdx = find(rem);
                if ~isempty(remIdx)
                    [~,o0] = sort(abs(sv(remIdx)),'ascend');
                    stab = remIdx(o0(1:min(k,numel(o0)))); % pct closest to zero
                else
                    stab = [];
                end

                lbl = nan(n,1); lbl(dec)=-1; lbl(inc)=1; lbl(stab)=0;
                labels(r,v)=lbl;
        end
    end
end

idx.decrease = labels==-1; idx.stable = labels==0; idx.increase = labels==1;

if flipBack
    labels=labels.'; fn=fieldnames(idx);
    for i=1:numel(fn), idx.(fn{i})=idx.(fn{i}).'; end
end
end

