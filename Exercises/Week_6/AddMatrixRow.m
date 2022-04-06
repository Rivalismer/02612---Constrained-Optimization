function A = AddMatrixRow(Ao, W)

A = [];
sw = length(W);
for i = 1:sw
    if W(i)
        A = [A Ao(:,i)];
    end
end
