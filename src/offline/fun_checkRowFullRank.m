function Is_RowFullRank = fun_checkRowFullRank(matrixCell)
    
    Is_RowFullRank = true;

    for i = 1:size(matrixCell,1)
        for j = 1:size(matrixCell,2)
            A = matrixCell{i,j};
            m = size(A,1);
            currentRank = rank(A);
            
            % 检查是否是行满秩
            if currentRank < m
                Is_RowFullRank = false;  % 发现一个不满秩就立即返回
                disp(currentRank);
                return;
            end
        end
    end

end