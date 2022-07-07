function globeff = compute_con_matrix(tot_commit2)
	data = readNPY(tot_commit2);
	data_norm = weight_conversion(data,'normalize');
	globeff = efficiency_wei(data_norm, 0);
end

