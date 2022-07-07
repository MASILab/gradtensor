function modularity = compute_con_matrix(tot_commit2)
	data = readNPY(tot_commit2);
	[~, modularity] = modularity_und(data,1);
end

