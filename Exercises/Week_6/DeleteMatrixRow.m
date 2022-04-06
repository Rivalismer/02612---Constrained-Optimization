function A = DeleteMatrixRow(Ao,num)

m = size(Ao(:,1));
T1 = Ao(:,1:(num - 1));
T2 = Ao(:,(num + 1):m(1));
A = [T1 T2];