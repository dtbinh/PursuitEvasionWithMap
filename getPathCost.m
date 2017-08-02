function Cost = getPathCost(Path)

[nx,~] = size(Path);

Cost = 0;
for k = 2:nx
    Cost = Cost + norm(Path(k,:) - Path(k-1,:));
end