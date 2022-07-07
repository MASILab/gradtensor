function charpathlen = compute_con_matrix(clen)
	data_len = readNPY(clen);
	distmatrix = distance_wei(data_len);
	[charpathlen,~, ~, ~, ~] = charpath(distmatrix, 1, 0);
end

