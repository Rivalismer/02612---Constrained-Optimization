function [x,lambda] = EqualityQPSolver(H,g,A,b,solver)

switch solver
    case 'LUSparse'
        [x,lambda] = EqualityQPSolverLUSparse(H,g,A,b);
   
    case 'LUDense'
        [x,lambda] = EqualityQPSolverLUDense(H,g,A,b);

    case 'LDLSparse'
        [x,lambda] = EqualityQPSolverLDLSparse(H,g,A,b);

    case 'LDLDense'
        [x,lambda] = EqualityQPSolverLDLDense(H,g,A,b);

    case 'NullSpace'
        [x,lambda] = EqualityQPSolverNullSpace(H,g,A,b);

    case 'RangeSpace'
        [x,lambda] = EqualityQPSolverRangeSpace(H,g,A,b);

    otherwise
        fprintf("Unknown solver");
        x = nan;
        lambda = nan;
end