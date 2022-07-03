function modularity,globeff,charpathlen = compute_con_matrix(tot_commit2,clen)
	data = readNPY(tot_commit2);
	[~, modularity] = modularity_und(data,1);
	%mod = [mod; modularity]

	data_norm = weight_conversion(data,'normalize');
	globeff = efficiency_wei(data_norm, 0);
	%ge = [ge; globeff];

	data_len = readNPY(clen);
	distmatrix = distance_wei(data_len);
	[charpathlen(i),~, ~, ~, ~] = charpath(distmatrix, 1, 0);
	%cpl = [cpl; charpathlen];
end

