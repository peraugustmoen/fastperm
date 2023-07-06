 #include "header.h"

 SEXP C_test_sparse_functions_with_example(){

	int n = 1000;


	
	SEXP alphaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP betaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP gammaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP logfactorialsSEXP = PROTECT(allocVector(REALSXP, n+1));
	SEXP mSEXP = PROTECT(allocVector(INTSXP, 1));
	SEXP k_SEXP = PROTECT(allocVector(INTSXP, 1));

	
	int * alpha = INTEGER(alphaSEXP);
	int * beta = INTEGER(betaSEXP);
	int * gamma = INTEGER(gammaSEXP);
	double * log_factorials = REAL(logfactorialsSEXP);
	int * m = INTEGER(mSEXP);
	int * k = INTEGER(k_SEXP);


	dictionary * new_log_subperms = init_dictionary(n);
	dictionary * old_log_subperms = init_dictionary(n);



	
	memset(alpha, 0, sizeof(int)*n);
	memset(beta, 0, sizeof(int)*n);
	memset(gamma, 0, sizeof(int)*n);
	memset(log_factorials, 0, sizeof(double)*(n+1));
	memset(m, 0, sizeof(int));
	memset(k, 0, sizeof(int));

	log_factorials[0]=0.0;
	//Rprintf("log_factorials[0] = %f\n", log_factorials[0]);
	for (int i = 1; i <= n; ++i)
	{
		log_factorials[i] = log_factorials[i-1] +log((double)(i));
		//Rprintf("log_factorials[%d] = %f\n", i, log_factorials[i]);
	}


	/*alpha[0] = 3634;
	alpha[1] = 17;
	alpha[2] = 82;
	alpha[3] = 14;
	alpha[4] = 2;
	alpha[5] = 2;
	alpha[6] = 4;
	alpha[7] = 1;
	alpha[8] = 1;
	alpha[9] = 38;
	alpha[10] = 399;
	alpha[11] = 88;
	alpha[12] = 3;
	alpha[13] = 105;
	alpha[14] = 2;
	alpha[15] = 2;
	alpha[16] = 811;
	alpha[17] = 238;
	alpha[18] = 10;
	alpha[19] = 3;
	alpha[20] = 4544;

	beta[0] = 0;
	beta[1] = 0;
	beta[2] = 9;
	beta[3] = 29;
	beta[4] = 13;
	beta[5] = 7;
	beta[6] = 7;
	beta[7] = 67;
	beta[8] = 144;
	beta[9] = 231;
	beta[10] = 64;
	beta[11] = 228;
	beta[12] = 479;
	beta[13] = 190;
	beta[14] = 212;
	beta[15] = 381;
	beta[16] = 761;
	beta[17] = 81;
	beta[18] = 472;
	beta[19] = 1579;

	gamma[0] = 2;
	gamma[1] = 50;
	gamma[2] = 251;
	gamma[3] = 586;
	gamma[4] = 109;
	gamma[5] = 84;
	gamma[6] = 96;
	gamma[7] = 410;
	gamma[8] = 562;
	gamma[9] = 548;
	gamma[10] = 106;
	gamma[11] = 406;
	gamma[12] = 598;
	gamma[13] = 172;
	gamma[14] = 176;
	gamma[15] = 240;
	gamma[16] = 336;
	gamma[17] = 25;
	gamma[18] = 135;
	gamma[19] = 154;

	*m = n;
	*k = 21;*/


	



	/*
	// n = 1000
	alpha[0] = 424;
	alpha[1] = 76;
	alpha[2] = 2;
	alpha[3] = 24;
	alpha[4] = 47;
	alpha[5] = 77;
	alpha[6] = 39;
	alpha[7] = 2;
	alpha[8] = 1;
	alpha[9] = 5;
	alpha[10] = 14;
	alpha[11] = 24;
	alpha[12] = 3;
	alpha[13] = 2;
	alpha[14] = 46;
	alpha[15] = 4;
	alpha[16] = 10;
	alpha[17] = 200;

	beta[0] = 0;
	beta[1] = 2;
	beta[2] = 5;
	beta[3] = 3;
	beta[4] = 28;
	beta[5] = 14;
	beta[6] = 5;
	beta[7] = 107;
	beta[8] = 22;
	beta[9] = 11;
	beta[10] = 9;
	beta[11] = 11;
	beta[12] = 56;
	beta[13] = 153;
	beta[14] = 72;
	beta[15] = 17;
	beta[16] = 5;

	gamma[0] = 44;
	gamma[1] = 39;
	gamma[2] = 29;
	gamma[3] = 12;
	gamma[4] = 88;
	gamma[5] = 36;
	gamma[6] = 8;
	gamma[7] = 127;
	gamma[8] = 23;
	gamma[9] = 13;
	gamma[10] = 7;
	gamma[11] = 6;
	gamma[12] = 13;
	gamma[13] = 30;
	gamma[14] = 5;
	gamma[15] = 0;
	gamma[16] = 0;

	*m = n;
	*k = 18;*/

alpha[0] = 3;
alpha[1] = 2;
alpha[2] = 1;
alpha[3] = 1;
alpha[4] = 1;
alpha[5] = 1;
alpha[6] = 1;
alpha[7] = 1;
alpha[8] = 2;
alpha[9] = 1;
alpha[10] = 2;
alpha[11] = 1;
alpha[12] = 1;
alpha[13] = 1;
alpha[14] = 1;
alpha[15] = 2;
alpha[16] = 4;
alpha[17] = 1;
alpha[18] = 1;
alpha[19] = 1;
alpha[20] = 1;
alpha[21] = 2;
alpha[22] = 1;
alpha[23] = 2;
alpha[24] = 1;
alpha[25] = 1;
alpha[26] = 2;
alpha[27] = 1;
alpha[28] = 2;
alpha[29] = 1;
alpha[30] = 1;
alpha[31] = 1;
alpha[32] = 3;
alpha[33] = 3;
alpha[34] = 2;
alpha[35] = 2;
alpha[36] = 1;
alpha[37] = 2;
alpha[38] = 1;
alpha[39] = 1;
alpha[40] = 1;
alpha[41] = 1;
alpha[42] = 1;
alpha[43] = 2;
alpha[44] = 2;
alpha[45] = 2;
alpha[46] = 1;
alpha[47] = 1;
alpha[48] = 1;
alpha[49] = 1;
alpha[50] = 1;
alpha[51] = 2;
alpha[52] = 1;
alpha[53] = 4;
alpha[54] = 1;
alpha[55] = 1;
alpha[56] = 1;
alpha[57] = 1;
alpha[58] = 1;
alpha[59] = 3;
alpha[60] = 1;
alpha[61] = 1;
alpha[62] = 1;
alpha[63] = 1;
alpha[64] = 2;
alpha[65] = 2;
alpha[66] = 1;
alpha[67] = 1;
alpha[68] = 1;
alpha[69] = 1;
alpha[70] = 2;
alpha[71] = 1;
alpha[72] = 1;
alpha[73] = 1;
alpha[74] = 1;
alpha[75] = 1;
alpha[76] = 1;
alpha[77] = 2;
alpha[78] = 1;
alpha[79] = 1;
alpha[80] = 1;
alpha[81] = 1;
alpha[82] = 1;
alpha[83] = 1;
alpha[84] = 3;
alpha[85] = 2;
alpha[86] = 2;
alpha[87] = 2;
alpha[88] = 2;
alpha[89] = 2;
alpha[90] = 4;
alpha[91] = 2;
alpha[92] = 1;
alpha[93] = 1;
alpha[94] = 1;
alpha[95] = 2;
alpha[96] = 2;
alpha[97] = 2;
alpha[98] = 1;
alpha[99] = 1;
alpha[100] = 1;
alpha[101] = 4;
alpha[102] = 2;
alpha[103] = 3;
alpha[104] = 1;
alpha[105] = 1;
alpha[106] = 2;
alpha[107] = 2;
alpha[108] = 1;
alpha[109] = 1;
alpha[110] = 1;
alpha[111] = 1;
alpha[112] = 2;
alpha[113] = 3;
alpha[114] = 1;
alpha[115] = 1;
alpha[116] = 1;
alpha[117] = 2;
alpha[118] = 3;
alpha[119] = 2;
alpha[120] = 1;
alpha[121] = 2;
alpha[122] = 3;
alpha[123] = 3;
alpha[124] = 1;
alpha[125] = 1;
alpha[126] = 1;
alpha[127] = 1;
alpha[128] = 1;
alpha[129] = 2;
alpha[130] = 1;
alpha[131] = 1;
alpha[132] = 1;
alpha[133] = 2;
alpha[134] = 3;
alpha[135] = 1;
alpha[136] = 1;
alpha[137] = 1;
alpha[138] = 1;
alpha[139] = 3;
alpha[140] = 1;
alpha[141] = 3;
alpha[142] = 1;
alpha[143] = 2;
alpha[144] = 1;
alpha[145] = 1;
alpha[146] = 1;
alpha[147] = 1;
alpha[148] = 1;
alpha[149] = 2;
alpha[150] = 3;
alpha[151] = 2;
alpha[152] = 1;
alpha[153] = 1;
alpha[154] = 1;
alpha[155] = 1;
alpha[156] = 1;
alpha[157] = 1;
alpha[158] = 1;
alpha[159] = 1;
alpha[160] = 3;
alpha[161] = 1;
alpha[162] = 1;
alpha[163] = 1;
alpha[164] = 2;
alpha[165] = 3;
alpha[166] = 2;
alpha[167] = 1;
alpha[168] = 1;
alpha[169] = 1;
alpha[170] = 2;
alpha[171] = 1;
alpha[172] = 1;
alpha[173] = 1;
alpha[174] = 2;
alpha[175] = 1;
alpha[176] = 3;
alpha[177] = 1;
alpha[178] = 1;
alpha[179] = 2;
alpha[180] = 1;
alpha[181] = 2;
alpha[182] = 2;
alpha[183] = 1;
alpha[184] = 1;
alpha[185] = 4;
alpha[186] = 1;
alpha[187] = 1;
alpha[188] = 1;
alpha[189] = 1;
alpha[190] = 1;
alpha[191] = 3;
alpha[192] = 1;
alpha[193] = 2;
alpha[194] = 1;
alpha[195] = 2;
alpha[196] = 2;
alpha[197] = 1;
alpha[198] = 1;
alpha[199] = 1;
alpha[200] = 1;
alpha[201] = 3;
alpha[202] = 1;
alpha[203] = 3;
alpha[204] = 1;
alpha[205] = 1;
alpha[206] = 4;
alpha[207] = 1;
alpha[208] = 1;
alpha[209] = 2;
alpha[210] = 1;
alpha[211] = 1;
alpha[212] = 3;
alpha[213] = 2;
alpha[214] = 1;
alpha[215] = 1;
alpha[216] = 3;
alpha[217] = 1;
alpha[218] = 1;
alpha[219] = 2;
alpha[220] = 1;
alpha[221] = 1;
alpha[222] = 1;
alpha[223] = 3;
alpha[224] = 1;
alpha[225] = 3;
alpha[226] = 1;
alpha[227] = 1;
alpha[228] = 2;
alpha[229] = 2;
alpha[230] = 2;
alpha[231] = 2;
alpha[232] = 1;
alpha[233] = 2;
alpha[234] = 2;
alpha[235] = 2;
alpha[236] = 1;
alpha[237] = 1;
alpha[238] = 2;
alpha[239] = 1;
alpha[240] = 1;
alpha[241] = 2;
alpha[242] = 2;
alpha[243] = 3;
alpha[244] = 1;
alpha[245] = 1;
alpha[246] = 3;
alpha[247] = 2;
alpha[248] = 1;
alpha[249] = 2;
alpha[250] = 1;
alpha[251] = 1;
alpha[252] = 1;
alpha[253] = 1;
alpha[254] = 1;
alpha[255] = 2;
alpha[256] = 1;
alpha[257] = 2;
alpha[258] = 2;
alpha[259] = 2;
alpha[260] = 1;
alpha[261] = 2;
alpha[262] = 2;
alpha[263] = 1;
alpha[264] = 1;
alpha[265] = 3;
alpha[266] = 1;
alpha[267] = 2;
alpha[268] = 1;
alpha[269] = 2;
alpha[270] = 1;
alpha[271] = 3;
alpha[272] = 2;
alpha[273] = 1;
alpha[274] = 2;
alpha[275] = 2;
alpha[276] = 1;
alpha[277] = 1;
alpha[278] = 2;
alpha[279] = 2;
alpha[280] = 3;
alpha[281] = 2;
alpha[282] = 5;
alpha[283] = 2;
alpha[284] = 2;
alpha[285] = 1;
alpha[286] = 3;
alpha[287] = 2;
alpha[288] = 1;
alpha[289] = 2;
alpha[290] = 3;
alpha[291] = 2;
alpha[292] = 2;
alpha[293] = 3;
alpha[294] = 2;
alpha[295] = 1;
alpha[296] = 1;
alpha[297] = 3;
alpha[298] = 2;
alpha[299] = 1;
alpha[300] = 1;
alpha[301] = 1;
alpha[302] = 1;
alpha[303] = 2;
alpha[304] = 1;
alpha[305] = 1;
alpha[306] = 1;
alpha[307] = 1;
alpha[308] = 2;
alpha[309] = 4;
alpha[310] = 1;
alpha[311] = 1;
alpha[312] = 1;
alpha[313] = 3;
alpha[314] = 2;
alpha[315] = 2;
alpha[316] = 1;
alpha[317] = 1;
alpha[318] = 3;
alpha[319] = 1;
alpha[320] = 1;
alpha[321] = 1;
alpha[322] = 1;
alpha[323] = 2;
alpha[324] = 1;
alpha[325] = 1;
alpha[326] = 1;
alpha[327] = 1;
alpha[328] = 4;
alpha[329] = 1;
alpha[330] = 1;
alpha[331] = 1;
alpha[332] = 1;
alpha[333] = 1;
alpha[334] = 3;
alpha[335] = 1;
alpha[336] = 1;
alpha[337] = 2;
alpha[338] = 2;
alpha[339] = 1;
alpha[340] = 1;
alpha[341] = 1;
alpha[342] = 3;
alpha[343] = 5;
alpha[344] = 2;
alpha[345] = 1;
alpha[346] = 2;
alpha[347] = 3;
alpha[348] = 1;
alpha[349] = 1;
alpha[350] = 1;
alpha[351] = 2;
alpha[352] = 3;
alpha[353] = 2;
alpha[354] = 2;
alpha[355] = 3;
alpha[356] = 1;
alpha[357] = 1;
alpha[358] = 1;
alpha[359] = 3;
alpha[360] = 2;
alpha[361] = 3;
alpha[362] = 1;
alpha[363] = 1;
alpha[364] = 2;
alpha[365] = 1;
alpha[366] = 2;
alpha[367] = 2;
alpha[368] = 1;
alpha[369] = 1;
alpha[370] = 1;
alpha[371] = 1;
alpha[372] = 2;
alpha[373] = 1;
alpha[374] = 4;
alpha[375] = 1;
alpha[376] = 1;
alpha[377] = 2;
alpha[378] = 1;
alpha[379] = 1;
alpha[380] = 2;
alpha[381] = 1;
alpha[382] = 1;
alpha[383] = 1;
alpha[384] = 1;
alpha[385] = 2;
alpha[386] = 2;
alpha[387] = 2;
alpha[388] = 1;
alpha[389] = 1;
alpha[390] = 1;
alpha[391] = 1;
alpha[392] = 1;
alpha[393] = 2;
alpha[394] = 1;
alpha[395] = 1;
alpha[396] = 2;
alpha[397] = 1;
alpha[398] = 3;
alpha[399] = 1;
alpha[400] = 1;
alpha[401] = 1;
alpha[402] = 4;
alpha[403] = 1;
alpha[404] = 2;
alpha[405] = 1;
alpha[406] = 2;
alpha[407] = 1;
alpha[408] = 3;
alpha[409] = 1;
alpha[410] = 1;
alpha[411] = 1;
alpha[412] = 3;
alpha[413] = 3;
alpha[414] = 1;
alpha[415] = 1;
alpha[416] = 3;
alpha[417] = 1;
alpha[418] = 1;
alpha[419] = 2;
alpha[420] = 1;
alpha[421] = 2;
alpha[422] = 1;
alpha[423] = 1;
alpha[424] = 2;
alpha[425] = 2;
alpha[426] = 1;
alpha[427] = 1;
alpha[428] = 1;
alpha[429] = 2;
alpha[430] = 1;
alpha[431] = 1;
alpha[432] = 3;
alpha[433] = 2;
alpha[434] = 2;
alpha[435] = 1;
alpha[436] = 2;
alpha[437] = 3;
alpha[438] = 1;
alpha[439] = 4;
alpha[440] = 2;
alpha[441] = 4;
alpha[442] = 1;
alpha[443] = 1;
alpha[444] = 1;
alpha[445] = 1;
alpha[446] = 1;
alpha[447] = 2;
alpha[448] = 2;
alpha[449] = 2;
alpha[450] = 3;
alpha[451] = 2;
alpha[452] = 2;
alpha[453] = 3;
alpha[454] = 3;
alpha[455] = 1;
alpha[456] = 1;
alpha[457] = 1;
alpha[458] = 3;
alpha[459] = 1;
alpha[460] = 2;
alpha[461] = 1;
alpha[462] = 1;
alpha[463] = 1;
alpha[464] = 1;
alpha[465] = 2;
alpha[466] = 1;
alpha[467] = 2;
alpha[468] = 1;
alpha[469] = 2;
alpha[470] = 1;
alpha[471] = 1;
alpha[472] = 1;
alpha[473] = 2;
alpha[474] = 1;
alpha[475] = 2;
alpha[476] = 2;
alpha[477] = 2;
alpha[478] = 2;
alpha[479] = 1;
alpha[480] = 1;
alpha[481] = 1;
alpha[482] = 2;
alpha[483] = 1;
alpha[484] = 1;
alpha[485] = 1;
alpha[486] = 1;
alpha[487] = 1;
alpha[488] = 1;
alpha[489] = 2;
alpha[490] = 1;
alpha[491] = 1;
alpha[492] = 2;
alpha[493] = 1;
alpha[494] = 1;
alpha[495] = 2;
alpha[496] = 2;
alpha[497] = 2;
alpha[498] = 3;
alpha[499] = 1;
alpha[500] = 2;
alpha[501] = 1;
alpha[502] = 1;
alpha[503] = 2;
alpha[504] = 2;
alpha[505] = 1;
alpha[506] = 1;
alpha[507] = 1;
alpha[508] = 1;
alpha[509] = 1;
alpha[510] = 2;
alpha[511] = 1;
alpha[512] = 1;
alpha[513] = 2;
alpha[514] = 1;
alpha[515] = 1;
alpha[516] = 1;
alpha[517] = 1;
alpha[518] = 1;
alpha[519] = 1;
alpha[520] = 4;
alpha[521] = 1;
alpha[522] = 1;
alpha[523] = 2;
alpha[524] = 3;
alpha[525] = 2;
alpha[526] = 2;
alpha[527] = 3;
alpha[528] = 3;
alpha[529] = 1;
alpha[530] = 2;
alpha[531] = 2;
alpha[532] = 3;
alpha[533] = 1;
alpha[534] = 1;
alpha[535] = 1;
alpha[536] = 2;
alpha[537] = 1;
alpha[538] = 1;
alpha[539] = 1;
alpha[540] = 1;
alpha[541] = 1;
alpha[542] = 2;
alpha[543] = 1;
alpha[544] = 1;
alpha[545] = 1;
alpha[546] = 1;
alpha[547] = 1;
alpha[548] = 1;
alpha[549] = 2;
alpha[550] = 1;
alpha[551] = 3;
alpha[552] = 1;
alpha[553] = 1;
alpha[554] = 1;
alpha[555] = 1;
alpha[556] = 1;
alpha[557] = 2;
alpha[558] = 2;
alpha[559] = 1;
alpha[560] = 1;
alpha[561] = 1;
alpha[562] = 1;
alpha[563] = 1;
alpha[564] = 3;
alpha[565] = 4;
alpha[566] = 1;
alpha[567] = 1;
alpha[568] = 1;
alpha[569] = 1;
alpha[570] = 2;
alpha[571] = 2;
alpha[572] = 1;
alpha[573] = 1;
alpha[574] = 1;
alpha[575] = 1;
alpha[576] = 2;
alpha[577] = 1;
alpha[578] = 3;
alpha[579] = 1;
alpha[580] = 1;
alpha[581] = 2;
alpha[582] = 1;
alpha[583] = 1;
alpha[584] = 2;
alpha[585] = 1;
alpha[586] = 1;
alpha[587] = 1;
alpha[588] = 1;
alpha[589] = 1;
alpha[590] = 3;
alpha[591] = 2;
alpha[592] = 1;
alpha[593] = 1;
alpha[594] = 2;
alpha[595] = 1;
alpha[596] = 2;
alpha[597] = 1;
alpha[598] = 3;
alpha[599] = 2;
alpha[600] = 1;
alpha[601] = 2;
alpha[602] = 1;
alpha[603] = 2;
alpha[604] = 1;
alpha[605] = 1;
alpha[606] = 1;
alpha[607] = 1;
alpha[608] = 4;
alpha[609] = 1;
alpha[610] = 1;
alpha[611] = 4;
alpha[612] = 1;
alpha[613] = 2;
alpha[614] = 1;
alpha[615] = 3;
alpha[616] = 1;
alpha[617] = 2;
alpha[618] = 1;
alpha[619] = 1;
alpha[620] = 3;
alpha[621] = 1;
alpha[622] = 1;
alpha[623] = 3;
alpha[624] = 3;
alpha[625] = 2;
alpha[626] = 1;
alpha[627] = 2;
alpha[628] = 1;
alpha[629] = 2;
alpha[630] = 2;
alpha[631] = 1;
alpha[632] = 1;
alpha[633] = 1;
alpha[634] = 2;

beta[0] = 0;
beta[1] = 0;
beta[2] = 0;
beta[3] = 0;
beta[4] = 0;
beta[5] = 0;
beta[6] = 0;
beta[7] = 0;
beta[8] = 0;
beta[9] = 0;
beta[10] = 0;
beta[11] = 0;
beta[12] = 0;
beta[13] = 0;
beta[14] = 0;
beta[15] = 0;
beta[16] = 0;
beta[17] = 0;
beta[18] = 0;
beta[19] = 0;
beta[20] = 0;
beta[21] = 0;
beta[22] = 0;
beta[23] = 0;
beta[24] = 0;
beta[25] = 0;
beta[26] = 0;
beta[27] = 0;
beta[28] = 0;
beta[29] = 0;
beta[30] = 0;
beta[31] = 0;
beta[32] = 0;
beta[33] = 0;
beta[34] = 0;
beta[35] = 0;
beta[36] = 0;
beta[37] = 0;
beta[38] = 0;
beta[39] = 0;
beta[40] = 0;
beta[41] = 0;
beta[42] = 0;
beta[43] = 0;
beta[44] = 0;
beta[45] = 0;
beta[46] = 0;
beta[47] = 0;
beta[48] = 0;
beta[49] = 0;
beta[50] = 0;
beta[51] = 0;
beta[52] = 0;
beta[53] = 0;
beta[54] = 0;
beta[55] = 0;
beta[56] = 0;
beta[57] = 0;
beta[58] = 0;
beta[59] = 0;
beta[60] = 0;
beta[61] = 0;
beta[62] = 0;
beta[63] = 0;
beta[64] = 0;
beta[65] = 0;
beta[66] = 0;
beta[67] = 0;
beta[68] = 0;
beta[69] = 0;
beta[70] = 0;
beta[71] = 0;
beta[72] = 0;
beta[73] = 0;
beta[74] = 0;
beta[75] = 0;
beta[76] = 0;
beta[77] = 0;
beta[78] = 0;
beta[79] = 0;
beta[80] = 0;
beta[81] = 0;
beta[82] = 0;
beta[83] = 0;
beta[84] = 0;
beta[85] = 0;
beta[86] = 0;
beta[87] = 0;
beta[88] = 0;
beta[89] = 0;
beta[90] = 0;
beta[91] = 0;
beta[92] = 0;
beta[93] = 0;
beta[94] = 0;
beta[95] = 0;
beta[96] = 0;
beta[97] = 1;
beta[98] = 0;
beta[99] = 0;
beta[100] = 0;
beta[101] = 0;
beta[102] = 0;
beta[103] = 0;
beta[104] = 0;
beta[105] = 0;
beta[106] = 0;
beta[107] = 0;
beta[108] = 0;
beta[109] = 0;
beta[110] = 0;
beta[111] = 0;
beta[112] = 0;
beta[113] = 1;
beta[114] = 0;
beta[115] = 0;
beta[116] = 0;
beta[117] = 0;
beta[118] = 0;
beta[119] = 0;
beta[120] = 0;
beta[121] = 3;
beta[122] = 0;
beta[123] = 0;
beta[124] = 0;
beta[125] = 0;
beta[126] = 0;
beta[127] = 0;
beta[128] = 0;
beta[129] = 1;
beta[130] = 0;
beta[131] = 0;
beta[132] = 0;
beta[133] = 1;
beta[134] = 0;
beta[135] = 0;
beta[136] = 0;
beta[137] = 0;
beta[138] = 0;
beta[139] = 0;
beta[140] = 1;
beta[141] = 1;
beta[142] = 0;
beta[143] = 1;
beta[144] = 1;
beta[145] = 1;
beta[146] = 1;
beta[147] = 1;
beta[148] = 0;
beta[149] = 0;
beta[150] = 0;
beta[151] = 0;
beta[152] = 1;
beta[153] = 0;
beta[154] = 0;
beta[155] = 0;
beta[156] = 0;
beta[157] = 0;
beta[158] = 0;
beta[159] = 0;
beta[160] = 0;
beta[161] = 0;
beta[162] = 0;
beta[163] = 0;
beta[164] = 0;
beta[165] = 1;
beta[166] = 0;
beta[167] = 0;
beta[168] = 0;
beta[169] = 0;
beta[170] = 0;
beta[171] = 0;
beta[172] = 0;
beta[173] = 0;
beta[174] = 0;
beta[175] = 0;
beta[176] = 0;
beta[177] = 0;
beta[178] = 0;
beta[179] = 0;
beta[180] = 2;
beta[181] = 0;
beta[182] = 0;
beta[183] = 1;
beta[184] = 0;
beta[185] = 0;
beta[186] = 0;
beta[187] = 0;
beta[188] = 0;
beta[189] = 0;
beta[190] = 0;
beta[191] = 1;
beta[192] = 0;
beta[193] = 0;
beta[194] = 1;
beta[195] = 0;
beta[196] = 0;
beta[197] = 0;
beta[198] = 0;
beta[199] = 0;
beta[200] = 1;
beta[201] = 0;
beta[202] = 1;
beta[203] = 0;
beta[204] = 0;
beta[205] = 1;
beta[206] = 0;
beta[207] = 1;
beta[208] = 1;
beta[209] = 0;
beta[210] = 0;
beta[211] = 1;
beta[212] = 0;
beta[213] = 0;
beta[214] = 1;
beta[215] = 0;
beta[216] = 0;
beta[217] = 0;
beta[218] = 0;
beta[219] = 0;
beta[220] = 0;
beta[221] = 1;
beta[222] = 0;
beta[223] = 0;
beta[224] = 0;
beta[225] = 0;
beta[226] = 0;
beta[227] = 0;
beta[228] = 0;
 
beta[229] = 0;
beta[230] = 1;
beta[231] = 0;
beta[232] = 0;
beta[233] = 0;
beta[234] = 1;
beta[235] = 0;
beta[236] = 0;
beta[237] = 0;
beta[238] = 1;
beta[239] = 0;
beta[240] = 1;
beta[241] = 0;
beta[242] = 1;
beta[243] = 0;
beta[244] = 0;
beta[245] = 1;
beta[246] = 0;
beta[247] = 1;
beta[248] = 1;
beta[249] = 0;
beta[250] = 0;
beta[251] = 0;
beta[252] = 0;
beta[253] = 0;
beta[254] = 0;
beta[255] = 0;
beta[256] = 1;
beta[257] = 0;
beta[258] = 0;
beta[259] = 0;
beta[260] = 0;
beta[261] = 1;
beta[262] = 0;
beta[263] = 1;
beta[264] = 0;
beta[265] = 0;
beta[266] = 0;
beta[267] = 0;
beta[268] = 1;
beta[269] = 0;
beta[270] = 2;
beta[271] = 0;
beta[272] = 3;
beta[273] = 0;
beta[274] = 1;
beta[275] = 1;
beta[276] = 4;
beta[277] = 0;
beta[278] = 0;
beta[279] = 0;
beta[280] = 1;
beta[281] = 1;
beta[282] = 1;
beta[283] = 1;
beta[284] = 2;
beta[285] = 1;
beta[286] = 0;
beta[287] = 1;
beta[288] = 0;
beta[289] = 2;
beta[290] = 0;
beta[291] = 0;
beta[292] = 1;
beta[293] = 1;
beta[294] = 0;
beta[295] = 0;
beta[296] = 1;
beta[297] = 0;
beta[298] = 0;
beta[299] = 1;
beta[300] = 0;
beta[301] = 1;
beta[302] = 1;
beta[303] = 1;
beta[304] = 2;
beta[305] = 0;
beta[306] = 0;
beta[307] = 0;
beta[308] = 0;
beta[309] = 1;
beta[310] = 1;
beta[311] = 0;
beta[312] = 1;
beta[313] = 1;
beta[314] = 0;
beta[315] = 1;
beta[316] = 0;
beta[317] = 0;
beta[318] = 0;
beta[319] = 1;
beta[320] = 1;
beta[321] = 1;
beta[322] = 0;
beta[323] = 0;
beta[324] = 1;
beta[325] = 2;
beta[326] = 1;
beta[327] = 0;
beta[328] = 1;
beta[329] = 0;
beta[330] = 1;
beta[331] = 1;
beta[332] = 0;
beta[333] = 1;
beta[334] = 1;
beta[335] = 1;
beta[336] = 2;
beta[337] = 0;
beta[338] = 2;
beta[339] = 0;
beta[340] = 2;
beta[341] = 0;
beta[342] = 1;
beta[343] = 0;
beta[344] = 1;
beta[345] = 1;
beta[346] = 0;
beta[347] = 0;
beta[348] = 1;
beta[349] = 0;
beta[350] = 1;
beta[351] = 2;
beta[352] = 3;
beta[353] = 1;
beta[354] = 1;
beta[355] = 1;
beta[356] = 1;
beta[357] = 1;
beta[358] = 3;
beta[359] = 1;
beta[360] = 2;
beta[361] = 1;
beta[362] = 0;
beta[363] = 1;
beta[364] = 1;
beta[365] = 2;
beta[366] = 4;
beta[367] = 1;
beta[368] = 1;
beta[369] = 2;
beta[370] = 3;
beta[371] = 1;
beta[372] = 3;
beta[373] = 1;
beta[374] = 0;
beta[375] = 2;
beta[376] = 1;
beta[377] = 1;
beta[378] = 0;
beta[379] = 1;
beta[380] = 2;
beta[381] = 2;
beta[382] = 1;
beta[383] = 0;
beta[384] = 2;
beta[385] = 0;
beta[386] = 1;
beta[387] = 0;
beta[388] = 1;
beta[389] = 1;
beta[390] = 0;
beta[391] = 0;
beta[392] = 1;
beta[393] = 3;
beta[394] = 2;
beta[395] = 1;
beta[396] = 1;
beta[397] = 1;
beta[398] = 2;
beta[399] = 1;
beta[400] = 1;
beta[401] = 1;
beta[402] = 1;
beta[403] = 1;
beta[404] = 1;
beta[405] = 1;
beta[406] = 1;
beta[407] = 1;
beta[408] = 0;
beta[409] = 3;
beta[410] = 1;
beta[411] = 1;
beta[412] = 1;
beta[413] = 0;
beta[414] = 0;
beta[415] = 1;
beta[416] = 1;
beta[417] = 1;
beta[418] = 3;
beta[419] = 2;
beta[420] = 1;
beta[421] = 4;
beta[422] = 1;
beta[423] = 1;
beta[424] = 0;
beta[425] = 1;
beta[426] = 2;
beta[427] = 1;
beta[428] = 2;
beta[429] = 3;
beta[430] = 1;
beta[431] = 1;
beta[432] = 2;
beta[433] = 0;
beta[434] = 1;
beta[435] = 1;
beta[436] = 2;
beta[437] = 1;
beta[438] = 2;
beta[439] = 1;
beta[440] = 1;
beta[441] = 1;
beta[442] = 1;
beta[443] = 2;
beta[444] = 1;
beta[445] = 2;
beta[446] = 0;
beta[447] = 2;
beta[448] = 1;
beta[449] = 1;
beta[450] = 0;
beta[451] = 2;
beta[452] = 1;
beta[453] = 1;
beta[454] = 1;
beta[455] = 1;
beta[456] = 1;
beta[457] = 1;
beta[458] = 1;
beta[459] = 3;
beta[460] = 1;
beta[461] = 1;
beta[462] = 1;
beta[463] = 1;
beta[464] = 1;
beta[465] = 0;
beta[466] = 2;
beta[467] = 1;
beta[468] = 1;
beta[469] = 2;
beta[470] = 1;
beta[471] = 1;
beta[472] = 3;
beta[473] = 3;
beta[474] = 1;
beta[475] = 2;
beta[476] = 1;
beta[477] = 2;
beta[478] = 2;
beta[479] = 1;
beta[480] = 3;
beta[481] = 0;
beta[482] = 1;
beta[483] = 2;
beta[484] = 2;
beta[485] = 1;
beta[486] = 1;
beta[487] = 2;
beta[488] = 1;
beta[489] = 2;
beta[490] = 1;
beta[491] = 1;
beta[492] = 4;
beta[493] = 6;
beta[494] = 1;
beta[495] = 1;
beta[496] = 1;
beta[497] = 1;
beta[498] = 1;
beta[499] = 1;
beta[500] = 1;
beta[501] = 1;
beta[502] = 2;
beta[503] = 1;
beta[504] = 2;
beta[505] = 1;
beta[506] = 2;
beta[507] = 1;
beta[508] = 2;
beta[509] = 3;
beta[510] = 1;
beta[511] = 3;
beta[512] = 3;
beta[513] = 1;
beta[514] = 1;
beta[515] = 1;
beta[516] = 1;
beta[517] = 1;
beta[518] = 2;
beta[519] = 1;
beta[520] = 1;
beta[521] = 2;
beta[522] = 1;
beta[523] = 4;
beta[524] = 1;
beta[525] = 1;
beta[526] = 1;
beta[527] = 3;
beta[528] = 2;
beta[529] = 2;
beta[530] = 1;
beta[531] = 1;
beta[532] = 1;
beta[533] = 3;
beta[534] = 2;
beta[535] = 1;
beta[536] = 1;
beta[537] = 1;
beta[538] = 1;
beta[539] = 4;
beta[540] = 3;
beta[541] = 1;
beta[542] = 1;
beta[543] = 4;
beta[544] = 1;
beta[545] = 2;
beta[546] = 1;
beta[547] = 0;
beta[548] = 1;
beta[549] = 1;
beta[550] = 1;
beta[551] = 4;
beta[552] = 3;
beta[553] = 1;
beta[554] = 1;
beta[555] = 1;
beta[556] = 1;
beta[557] = 2;
beta[558] = 3;
beta[559] = 1;
beta[560] = 1;
beta[561] = 1;
beta[562] = 1;
beta[563] = 1;
beta[564] = 1;
beta[565] = 1;
beta[566] = 2;
beta[567] = 2;
beta[568] = 1;
beta[569] = 3;
beta[570] = 1;
beta[571] = 2;
beta[572] = 1;
beta[573] = 1;
beta[574] = 2;
beta[575] = 4;
beta[576] = 1;
beta[577] = 2;
beta[578] = 3;
beta[579] = 1;
beta[580] = 4;
beta[581] = 1;
beta[582] = 1;
beta[583] = 4;
beta[584] = 1;
beta[585] = 1;
beta[586] = 1;
beta[587] = 1;
beta[588] = 1;
beta[589] = 1;
beta[590] = 3;
beta[591] = 1;
beta[592] = 1;
beta[593] = 1;
beta[594] = 1;
beta[595] = 1;
beta[596] = 1;
beta[597] = 1;
beta[598] = 1;
beta[599] = 1;
beta[600] = 1;
beta[601] = 1;
beta[602] = 1;
beta[603] = 2;
beta[604] = 1;
beta[605] = 1;
beta[606] = 1;
beta[607] = 1;
beta[608] = 3;
beta[609] = 2;
beta[610] = 2;
beta[611] = 1;
beta[612] = 1;
beta[613] = 1;
beta[614] = 1;
beta[615] = 2;
beta[616] = 1;
beta[617] = 2;
beta[618] = 1;
beta[619] = 1;
beta[620] = 1;
beta[621] = 1;
beta[622] = 6;
beta[623] = 4;
beta[624] = 1;
beta[625] = 1;
beta[626] = 2;
beta[627] = 2;
beta[628] = 2;
beta[629] = 2;
beta[630] = 1;
beta[631] = 1;
beta[632] = 2;
beta[633] = 1;

gamma[0] = 2;
gamma[1] = 3;
gamma[2] = 1;
gamma[3] = 2;
gamma[4] = 1;
gamma[5] = 1;
gamma[6] = 3;
gamma[7] = 1;
gamma[8] = 1;
gamma[9] = 1;
gamma[10] = 1;
gamma[11] = 1;
gamma[12] = 2;
gamma[13] = 3;
gamma[14] = 1;
gamma[15] = 1;
gamma[16] = 2;
gamma[17] = 1;
gamma[18] = 1;
gamma[19] = 2;
gamma[20] = 2;
gamma[21] = 1;
gamma[22] = 1;
gamma[23] = 1;
gamma[24] = 1;
gamma[25] = 1;
gamma[26] = 3;
gamma[27] = 1;
gamma[28] = 1;
gamma[29] = 1;
gamma[30] = 1;
gamma[31] = 2;
gamma[32] = 1;
gamma[33] = 1;
gamma[34] = 1;
gamma[35] = 3;
gamma[36] = 4;
gamma[37] = 1;
gamma[38] = 2;
gamma[39] = 1;
gamma[40] = 1;
gamma[41] = 2;
gamma[42] = 1;
gamma[43] = 1;
gamma[44] = 1;
gamma[45] = 1;
gamma[46] = 1;
gamma[47] = 1;
gamma[48] = 2;
gamma[49] = 1;
gamma[50] = 1;
gamma[51] = 2;
gamma[52] = 3;
gamma[53] = 2;
gamma[54] = 1;
gamma[55] = 1;
gamma[56] = 1;
gamma[57] = 1;
gamma[58] = 1;
gamma[59] = 1;
gamma[60] = 1;
gamma[61] = 1;
gamma[62] = 1;
gamma[63] = 1;
gamma[64] = 1;
gamma[65] = 2;
gamma[66] = 1;
gamma[67] = 1;
gamma[68] = 4;
gamma[69] = 1;
gamma[70] = 1;
gamma[71] = 1;
gamma[72] = 1;
gamma[73] = 1;
gamma[74] = 1;
gamma[75] = 2;
gamma[76] = 1;
gamma[77] = 1;
gamma[78] = 1;
gamma[79] = 1;
gamma[80] = 2;
gamma[81] = 1;
gamma[82] = 6;
gamma[83] = 1;
gamma[84] = 3;
gamma[85] = 1;
gamma[86] = 1;
gamma[87] = 1;
gamma[88] = 1;
gamma[89] = 1;
gamma[90] = 2;
gamma[91] = 1;
gamma[92] = 1;
gamma[93] = 1;
gamma[94] = 1;
gamma[95] = 1;
gamma[96] = 1;
gamma[97] = 2;
gamma[98] = 1;
gamma[99] = 1;
gamma[100] = 1;
gamma[101] = 1;
gamma[102] = 1;
gamma[103] = 2;
gamma[104] = 1;
gamma[105] = 1;
gamma[106] = 3;
gamma[107] = 1;
gamma[108] = 1;
gamma[109] = 1;
gamma[110] = 1;
gamma[111] = 2;
gamma[112] = 1;
gamma[113] = 2;
gamma[114] = 3;
gamma[115] = 1;
gamma[116] = 1;
gamma[117] = 1;
gamma[118] = 1;
gamma[119] = 1;
gamma[120] = 1;
gamma[121] = 0;
gamma[122] = 2;
gamma[123] = 1;
gamma[124] = 1;
gamma[125] = 3;
gamma[126] = 1;
gamma[127] = 2;
gamma[128] = 1;
gamma[129] = 1;
gamma[130] = 1;
gamma[131] = 2;
gamma[132] = 1;
gamma[133] = 2;
gamma[134] = 1;
gamma[135] = 1;
gamma[136] = 1;
gamma[137] = 1;
gamma[138] = 2;
gamma[139] = 1;
gamma[140] = 4;
gamma[141] = 1;
gamma[142] = 1;
gamma[143] = 0;
gamma[144] = 0;
gamma[145] = 0;
gamma[146] = 0;
gamma[147] = 1;
gamma[148] = 1;
gamma[149] = 1;
gamma[150] = 4;
gamma[151] = 2;
gamma[152] = 2;
gamma[153] = 1;
gamma[154] = 2;
gamma[155] = 4;
gamma[156] = 2;
gamma[157] = 1;
gamma[158] = 1;
gamma[159] = 1;
gamma[160] = 1;
gamma[161] = 1;
gamma[162] = 1;
gamma[163] = 1;
gamma[164] = 4;
gamma[165] = 1;
gamma[166] = 1;
gamma[167] = 1;
gamma[168] = 1;
gamma[169] = 1;
gamma[170] = 1;
gamma[171] = 3;
gamma[172] = 2;
gamma[173] = 1;
gamma[174] = 1;
gamma[175] = 1;
gamma[176] = 1;
gamma[177] = 4;
gamma[178] = 1;
gamma[179] = 2;
gamma[180] = 1;
gamma[181] = 4;
gamma[182] = 2;
gamma[183] = 1;
gamma[184] = 1;
gamma[185] = 1;
gamma[186] = 1;
gamma[187] = 1;
gamma[188] = 1;
gamma[189] = 2;
gamma[190] = 1;
gamma[191] = 0;
gamma[192] = 1;
gamma[193] = 1;
gamma[194] = 1;
gamma[195] = 1;
gamma[196] = 1;
gamma[197] = 1;
gamma[198] = 3;
gamma[199] = 2;
gamma[200] = 4;
gamma[201] = 1;
gamma[202] = 3;
gamma[203] = 1;
gamma[204] = 1;
gamma[205] = 0;
gamma[206] = 1;
gamma[207] = 2;
gamma[208] = 0;
gamma[209] = 1;
gamma[210] = 1;
gamma[211] = 0;
gamma[212] = 1;
gamma[213] = 1;
gamma[214] = 1;
gamma[215] = 2;
gamma[216] = 1;
gamma[217] = 3;
gamma[218] = 3;
gamma[219] = 1;
gamma[220] = 1;
gamma[221] = 2;
gamma[222] = 3;
gamma[223] = 1;
gamma[224] = 1;
gamma[225] = 1;
gamma[226] = 1;
gamma[227] = 1;
gamma[228] = 1;
gamma[229] = 2;
gamma[230] = 0;
gamma[231] = 1;
gamma[232] = 1;
gamma[233] = 1;
gamma[234] = 0;
gamma[235] = 1;
gamma[236] = 1;
gamma[237] = 2;
gamma[238] = 0;
gamma[239] = 3;
gamma[240] = 0;
gamma[241] = 2;
gamma[242] = 1;
gamma[243] = 1;
gamma[244] = 1;
gamma[245] = 0;
gamma[246] = 3;
gamma[247] = 0;
gamma[248] = 0;
gamma[249] = 1;
gamma[250] = 1;
gamma[251] = 1;
gamma[252] = 1;
gamma[253] = 3;
gamma[254] = 4;
gamma[255] = 1;
gamma[256] = 1;
gamma[257] = 2;
gamma[258] = 2;
gamma[259] = 1;
gamma[260] = 1;
gamma[261] = 0;
gamma[262] = 1;
gamma[263] = 0;
gamma[264] = 1;
gamma[265] = 1;
gamma[266] = 1;
gamma[267] = 1;
gamma[268] = 0;
gamma[269] = 1;
gamma[270] = 1;
gamma[271] = 2;
gamma[272] = 0;
gamma[273] = 1;
gamma[274] = 0;
gamma[275] = 0;
gamma[276] = 1;
gamma[277] = 1;
gamma[278] = 1;
gamma[279] = 1;
gamma[280] = 1;
gamma[281] = 0;
gamma[282] = 2;
gamma[283] = 0;
gamma[284] = 0;
gamma[285] = 0;
gamma[286] = 1;
gamma[287] = 1;
gamma[288] = 1;
gamma[289] = 1;
gamma[290] = 1;
gamma[291] = 1;
gamma[292] = 3;
gamma[293] = 0;
gamma[294] = 1;
gamma[295] = 1;
gamma[296] = 0;
gamma[297] = 1;
gamma[298] = 1;
gamma[299] = 1;
gamma[300] = 1;
gamma[301] = 0;
gamma[302] = 0;
gamma[303] = 0;
gamma[304] = 1;
gamma[305] = 1;
gamma[306] = 3;
gamma[307] = 1;
gamma[308] = 1;
gamma[309] = 0;
gamma[310] = 1;
gamma[311] = 6;
gamma[312] = 0;
gamma[313] = 0;
gamma[314] = 1;
gamma[315] = 0;
gamma[316] = 1;
gamma[317] = 2;
gamma[318] = 1;
gamma[319] = 0;
gamma[320] = 0;
gamma[321] = 1;
gamma[322] = 3;
gamma[323] = 1;
gamma[324] = 0;
gamma[325] = 1;
gamma[326] = 0;
gamma[327] = 1;
gamma[328] = 1;
gamma[329] = 1;
gamma[330] = 0;
gamma[331] = 1;
gamma[332] = 1;
gamma[333] = 0;
gamma[334] = 0;
gamma[335] = 2;
gamma[336] = 1;
gamma[337] = 1;
gamma[338] = 1;
gamma[339] = 1;
gamma[340] = 0;
gamma[341] = 1;
gamma[342] = 0;
gamma[343] = 1;
gamma[344] = 0;
gamma[345] = 1;
gamma[346] = 2;
gamma[347] = 1;
gamma[348] = 0;
gamma[349] = 1;
gamma[350] = 0;
gamma[351] = 1;
gamma[352] = 1;
gamma[353] = 0;
gamma[354] = 1;
gamma[355] = 0;
gamma[356] = 0;
gamma[357] = 0;
gamma[358] = 0;
gamma[359] = 1;
gamma[360] = 0;
gamma[361] = 0;
gamma[362] = 1;
gamma[363] = 0;
gamma[364] = 0;
gamma[365] = 1;
gamma[366] = 1;
gamma[367] = 0;
gamma[368] = 1;
gamma[369] = 0;
gamma[370] = 0;
gamma[371] = 0;
gamma[372] = 0;
gamma[373] = 0;
gamma[374] = 1;
gamma[375] = 0;
gamma[376] = 0;
gamma[377] = 1;
gamma[378] = 1;
gamma[379] = 1;
gamma[380] = 0;
gamma[381] = 0;
gamma[382] = 0;
gamma[383] = 2;
gamma[384] = 0;
gamma[385] = 1;
gamma[386] = 0;
gamma[387] = 1;
gamma[388] = 0;
gamma[389] = 0;
gamma[390] = 1;
gamma[391] = 1;
gamma[392] = 0;
gamma[393] = 1;
gamma[394] = 0;
gamma[395] = 0;
gamma[396] = 0;
gamma[397] = 0;
gamma[398] = 0;
gamma[399] = 0;
gamma[400] = 0;
gamma[401] = 1;
gamma[402] = 0;
gamma[403] = 0;
gamma[404] = 0;
gamma[405] = 0;
gamma[406] = 0;
gamma[407] = 0;
gamma[408] = 2;
gamma[409] = 0;
gamma[410] = 0;
gamma[411] = 0;
gamma[412] = 0;
gamma[413] = 1;
gamma[414] = 1;
gamma[415] = 0;
gamma[416] = 0;
gamma[417] = 2;
gamma[418] = 0;
gamma[419] = 0;
gamma[420] = 0;
gamma[421] = 0;
gamma[422] = 0;
gamma[423] = 0;
gamma[424] = 1;
gamma[425] = 0;
gamma[426] = 0;
gamma[427] = 0;
gamma[428] = 0;
gamma[429] = 1;
gamma[430] = 1;
gamma[431] = 0;
gamma[432] = 0;
gamma[433] = 1;
gamma[434] = 0;
gamma[435] = 0;
gamma[436] = 1;
gamma[437] = 0;
gamma[438] = 0;
gamma[439] = 0;
gamma[440] = 0;
gamma[441] = 0;
gamma[442] = 0;
gamma[443] = 0;
gamma[444] = 0;
gamma[445] = 1;
gamma[446] = 1;
gamma[447] = 0;
gamma[448] = 0;
gamma[449] = 0;
gamma[450] = 1;
gamma[451] = 0;
gamma[452] = 0;
gamma[453] = 1;
gamma[454] = 0;
gamma[455] = 0;
gamma[456] = 0;
gamma[457] = 0;
gamma[458] = 0;
gamma[459] = 0;
gamma[460] = 0;
gamma[461] = 1;
gamma[462] = 0;
gamma[463] = 0;
gamma[464] = 0;
gamma[465] = 1;
gamma[466] = 0;
gamma[467] = 0;
gamma[468] = 0;
gamma[469] = 0;
gamma[470] = 1;
gamma[471] = 1;
gamma[472] = 0;
gamma[473] = 0;
gamma[474] = 1;
gamma[475] = 0;
gamma[476] = 0;
gamma[477] = 0;
gamma[478] = 0;
gamma[479] = 0;
gamma[480] = 1;
gamma[481] = 1;
gamma[482] = 0;
gamma[483] = 0;
gamma[484] = 0;
gamma[485] = 0;
gamma[486] = 0;
gamma[487] = 0;
gamma[488] = 0;
gamma[489] = 1;
gamma[490] = 0;
gamma[491] = 0;
gamma[492] = 0;
gamma[493] = 1;
gamma[494] = 0;
gamma[495] = 0;
gamma[496] = 0;
gamma[497] = 1;
gamma[498] = 0;
gamma[499] = 0;
gamma[500] = 0;
gamma[501] = 0;
gamma[502] = 1;
gamma[503] = 0;
gamma[504] = 0;
gamma[505] = 0;
gamma[506] = 0;
gamma[507] = 0;
gamma[508] = 0;
gamma[509] = 0;
gamma[510] = 0;
gamma[511] = 1;
gamma[512] = 1;
gamma[513] = 0;
gamma[514] = 0;
gamma[515] = 0;
gamma[516] = 0;
gamma[517] = 0;
gamma[518] = 0;
gamma[519] = 0;
gamma[520] = 0;
gamma[521] = 0;
gamma[522] = 0;
gamma[523] = 0;
gamma[524] = 0;
gamma[525] = 0;
gamma[526] = 0;
gamma[527] = 0;
gamma[528] = 0;
gamma[529] = 0;
gamma[530] = 0;
gamma[531] = 0;
gamma[532] = 0;
gamma[533] = 0;
gamma[534] = 0;
gamma[535] = 0;
gamma[536] = 0;
gamma[537] = 0;
gamma[538] = 0;
gamma[539] = 0;
gamma[540] = 0;
gamma[541] = 0;
gamma[542] = 0;
gamma[543] = 2;
gamma[544] = 0;
gamma[545] = 0;
gamma[546] = 0;
gamma[547] = 1;
gamma[548] = 0;
gamma[549] = 0;
gamma[550] = 1;
gamma[551] = 0;
gamma[552] = 0;
gamma[553] = 0;
gamma[554] = 0;
gamma[555] = 0;
gamma[556] = 0;
gamma[557] = 0;
gamma[558] = 0;
gamma[559] = 0;
gamma[560] = 0;
gamma[561] = 0;
gamma[562] = 0;
gamma[563] = 0;
gamma[564] = 0;
gamma[565] = 0;
gamma[566] = 0;
gamma[567] = 0;
gamma[568] = 0;
gamma[569] = 0;
gamma[570] = 0;
gamma[571] = 0;
gamma[572] = 0;
gamma[573] = 0;
gamma[574] = 0;
gamma[575] = 0;
gamma[576] = 0;
gamma[577] = 0;
gamma[578] = 0;
gamma[579] = 0;
gamma[580] = 0;
gamma[581] = 0;
gamma[582] = 0;
gamma[583] = 0;
gamma[584] = 0;
gamma[585] = 0;
gamma[586] = 0;
gamma[587] = 0;
gamma[588] = 0;
gamma[589] = 0;
gamma[590] = 0;
gamma[591] = 0;
gamma[592] = 0;
gamma[593] = 0;
gamma[594] = 0;
gamma[595] = 0;
gamma[596] = 0;
gamma[597] = 0;
gamma[598] = 0;
gamma[599] = 0;
gamma[600] = 0;
gamma[601] = 0;
gamma[602] = 0;
gamma[603] = 0;
gamma[604] = 0;
gamma[605] = 0;
gamma[606] = 0;
gamma[607] = 0;
gamma[608] = 0;
gamma[609] = 0;
gamma[610] = 0;
gamma[611] = 0;
gamma[612] = 0;
gamma[613] = 0;
gamma[614] = 0;
gamma[615] = 0;
gamma[616] = 0;
gamma[617] = 0;
gamma[618] = 0;
gamma[619] = 0;
gamma[620] = 0;
gamma[621] = 0;
gamma[622] = 0;
gamma[623] = 0;
gamma[624] = 0;
gamma[625] = 0;
gamma[626] = 0;
gamma[627] = 0;
gamma[628] = 0;
gamma[629] = 0;
gamma[630] = 0;
gamma[631] = 0;
gamma[632] = 0;
gamma[633] = 0;

*m = 1000;
*k = 635;


	int history_len = 0;

	SEXP historySEXP = PROTECT(allocVector(INTSXP, 3*n));
	SEXP amount_historySEXP = PROTECT(allocVector(INTSXP, 6*n));

	int * history = INTEGER(historySEXP);
	int * amount_history = INTEGER(amount_historySEXP);

	memset(history, 0, sizeof(int)*3*n);
	memset(amount_history, 0, sizeof(int)*6*n);

/*	int t[2];
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			t[0] = i;
		    t[1] = j;
		    XXH64_hash_t hash = abs(XXH64(t, sizeof(int)*2, 1337));
		    Rprintf("hash: %d\n", ((int) hash) % n);
		}
	}
    
    return 0;*/

	Rprintf("REDUCING NOW\n");
	reduction(alpha,  beta,  gamma, m, n, k, history,
			   amount_history, &history_len, 0);

	Rprintf("history len = %d\n", history_len);

	Rprintf("REDUCED SUBPERMS\n");
	sparse_get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
					log_factorials, n,  m, k);
	//Rprintf("RESULT:\n");
	//print_matrix(n+1, n+1, new_log_subperms);

	dictionary * tmp  = old_log_subperms;
	old_log_subperms = new_log_subperms;
	new_log_subperms = tmp;


	/*UNPROTECT(11);
	return(mSEXP);*/

	Rprintf("==========\nReverse reduction:\n==========\n");
	//Rprintf("old = %d\n", old_log_subperms);
	//Rprintf("new = %d\n", new_log_subperms);
	dictionary * the_log_subperms = sparse_reverse_reduction(old_log_subperms, new_log_subperms, alpha,
					   beta,  gamma, m,  n, k,  history,
			           amount_history, &history_len, log_factorials);

	

	/*for (int i = 0; i < history_len; ++i)
	{
		Rprintf("amounthistory[%d] = (%d,%d)\n", i, amount_history[2*i], amount_history[2*i+1]);
	}*/

	
	double logperm =  Csparse_log_sum_exp(the_log_subperms);
	Rprintf("logperm = %f\n", logperm);



	UNPROTECT(8);
	free_dictionary(new_log_subperms);
	free_dictionary(old_log_subperms);




	return(mSEXP);

}



// THIS WORKS:
// SEXP C_test_functions_with_example(){
// 
// 	int n = 100;
// 
// 
// 	SEXP old_log_subpermsSEXP = PROTECT(allocVector(REALSXP, (n+1)*(n+1)));
// 	SEXP new_log_subpermsSEXP = PROTECT(allocVector(REALSXP, (n+1)*(n+1)));
// 	SEXP temp_vecSEXP = PROTECT(allocVector(REALSXP, n+1));
// 	SEXP alphaSEXP = PROTECT(allocVector(INTSXP, n));
// 	SEXP betaSEXP = PROTECT(allocVector(INTSXP, n));
// 	SEXP gammaSEXP = PROTECT(allocVector(INTSXP, n));
// 	SEXP logfactorialsSEXP = PROTECT(allocVector(REALSXP, n+1));
// 	SEXP mSEXP = PROTECT(allocVector(INTSXP, 1));
// 	SEXP k_SEXP = PROTECT(allocVector(INTSXP, 1));
// 
// 	double * old_log_subperms = REAL(old_log_subpermsSEXP);
// 	double * new_log_subperms = REAL(new_log_subpermsSEXP);
// 	double * temp_vec = REAL(temp_vecSEXP);
// 	int * alpha = INTEGER(alphaSEXP);
// 	int * beta = INTEGER(betaSEXP);
// 	int * gamma = INTEGER(gammaSEXP);
// 	double * log_factorials = REAL(logfactorialsSEXP);
// 	int * m = INTEGER(mSEXP);
// 	int * k = INTEGER(k_SEXP);
// 
// 
// 
// 
// 
// 
// 	memset(old_log_subperms, 0, sizeof(double)*(n+1)*(n+1));
// 	memset(old_log_subperms, 0, sizeof(double)*(n+1)*(n+1));
// 	memset(temp_vec, 0, sizeof(double)*(n+1));
// 	memset(alpha, 0, sizeof(int)*n);
// 	memset(beta, 0, sizeof(int)*n);
// 	memset(gamma, 0, sizeof(int)*n);
// 	memset(log_factorials, 0, sizeof(double)*(n+1));
// 	memset(m, 0, sizeof(int));
// 	memset(k, 0, sizeof(int));
// 
// 	log_factorials[0]=0.0;
// 	//Rprintf("log_factorials[0] = %f\n", log_factorials[0]);
// 	for (int i = 1; i <= n; ++i)
// 	{
// 		log_factorials[i] = log_factorials[i-1] +log((double)(i));
// 		//Rprintf("log_factorials[%d] = %f\n", i, log_factorials[i]);
// 	}
// 
// 
// 
// 
// 	for (int i = 0; i < n+1; ++i)
// 	{
// 		for (int j = 0; j < n+1; ++j)
// 		{
// 			new_log_subperms[cord_spec(i,j,(n+1))] = -1;
// 			old_log_subperms[cord_spec(i,j,(n+1))] = -1;
// 		}
// 	}
// 
// 
// 
// 	alpha[0] = 23;
// 	alpha[1] = 10;
// 	alpha[2] = 17;
// 	alpha[3] = 31;
// 	alpha[4] = 12;
// 	alpha[5] = 6;
// 	alpha[6] = 1;
// 
// 
// 	beta[0] = 2;
// 	beta[1] = 1;
// 	beta[2] = 4;
// 	beta[3] = 6;
// 	beta[4] = 9;
// 	beta[5] = 10;
// 
// 	gamma[0] =8;
// 	gamma[1] =9;
// 	gamma[2] =6;
// 	gamma[3] =4;
// 	gamma[4] =1;
// 	gamma[5] =0;
// 
// 	*m = 100;
// 	*k = 7;
// 
// 
// 	int history_len = 0;
// 
// 	SEXP historySEXP = PROTECT(allocVector(INTSXP, 3*n));
// 	SEXP amount_historySEXP = PROTECT(allocVector(INTSXP, 6*n));
// 
// 	int * history = INTEGER(historySEXP);
// 	int * amount_history = INTEGER(amount_historySEXP);
// 
// 	memset(history, 0, sizeof(int)*3*n);
// 	memset(amount_history, 0, sizeof(int)*6*n);
// 
// 	Rprintf("REDUCING NOW\n");
// 	reduction(alpha,  beta,  gamma, m, n, k, history,
// 			   amount_history, &history_len);
// 
// 
// 
// 	//Rprintf("REDUCED SUBPERMS\n");
// 	get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
// 					log_factorials, n,  m, k);
// 	//Rprintf("RESULT:\n");
// 	//print_matrix(n+1, n+1, new_log_subperms);
// 
// 	double * tmp  = old_log_subperms;
// 	old_log_subperms = new_log_subperms;
// 	new_log_subperms = tmp;
// 
// 
// 	/*UNPROTECT(11);
// 	return(mSEXP);*/
// 
// 	Rprintf("==========\nReverse reduction:\n==========\n");
// 	//Rprintf("old = %d\n", old_log_subperms);
// 	//Rprintf("new = %d\n", new_log_subperms);
// 	double * the_log_subperms = reverse_reduction(old_log_subperms, new_log_subperms, alpha,
// 					   beta,  gamma, m,  n, k,  history,
// 			           amount_history, &history_len, log_factorials,
// 			           temp_vec);
// 
// 	int maxxval = -1;
// 
// 	for (int i = 0; i < history_len; ++i)
// 	{
// 		Rprintf("amounthistory[%d] = (%d,%d)\n", i, amount_history[2*i], amount_history[2*i+1]);
// 	}
// 
// 	for (int i = 0; i < (n+1)*(n+1); ++i)
// 	{
// 		if(the_log_subperms[i]>maxxval){
// 			maxxval = the_log_subperms[i];
// 		}
// 	}
// 	double logperm =  Clog_sum_exp(the_log_subperms, (n+1)*(n+1),maxxval);
// 	Rprintf("logperm = %f\n", logperm);
// 
// 
// 
// 	UNPROTECT(11);
// 
// 
// 
// 	return(mSEXP);
// 
// }












// 
// SEXP C_test_sparse_functions_with_example(){
// 
// 	int n = 7;
// 
// 
// 	
// 	SEXP temp_vecSEXP = PROTECT(allocVector(REALSXP, n+1));
// 	SEXP alphaSEXP = PROTECT(allocVector(INTSXP, n));
// 	SEXP betaSEXP = PROTECT(allocVector(INTSXP, n));
// 	SEXP gammaSEXP = PROTECT(allocVector(INTSXP, n));
// 	SEXP logfactorialsSEXP = PROTECT(allocVector(REALSXP, n+1));
// 	SEXP mSEXP = PROTECT(allocVector(INTSXP, 1));
// 	SEXP k_SEXP = PROTECT(allocVector(INTSXP, 1));
// 
// 
// 	double * temp_vec = REAL(temp_vecSEXP);
// 	int * alpha = INTEGER(alphaSEXP);
// 	int * beta = INTEGER(betaSEXP);
// 	int * gamma = INTEGER(gammaSEXP);
// 	double * log_factorials = REAL(logfactorialsSEXP);
// 	int * m = INTEGER(mSEXP);
// 	int * k = INTEGER(k_SEXP);
// 
// 
// 
// 
// 
// 
// 	
// 	memset(temp_vec, 0, sizeof(double)*(n+1));
// 	memset(alpha, 0, sizeof(int)*n);
// 	memset(beta, 0, sizeof(int)*n);
// 	memset(gamma, 0, sizeof(int)*n);
// 	memset(log_factorials, 0, sizeof(double)*(n+1));
// 	memset(m, 0, sizeof(int));
// 	memset(k, 0, sizeof(int));
// 
// 	log_factorials[0]=0.0;
// 	Rprintf("log_factorials[0] = %f\n", log_factorials[0]);
// 	for (int i = 1; i <= n; ++i)
// 	{
// 		log_factorials[i] = log_factorials[i-1] +log((double)(i));
// 		Rprintf("log_factorials[%d] = %f\n", i, log_factorials[i]);
// 	}
// 
// 	dictionary * new_log_subperms = init_dictionary(10);
// 	dictionary * old_log_subperms = init_dictionary(10);
//   
//   
//   
// 	alpha[0] = 5;
// 	alpha[1] = 2;
// 
// 	beta[0] = 2;
// 	gamma[0] =0;
// 	(*m) = 2;
// 	*k = 2;
// 
// 	pair pp = {.x = 2, .y = 0};
// 	add_to_dictionary(old_log_subperms, pp, log(20));
// 	
// 
// 	Rprintf("Staring vals of alpha: ");
// 	print_int_vector(n, alpha);
// 
// 	Rprintf("after reverse bottom trim:\n");
// 	sparse_reverse_bt(old_log_subperms, new_log_subperms, 1,  0,
// 				alpha, beta,  gamma,  log_factorials, n, m,
// 				k);
// 	
// 	print_sparse_matrix(new_log_subperms);
// 
// 	Rprintf("Vectors:\n Alpha: ");
// 	print_int_vector(n, alpha);
// 	Rprintf("Beta:");
// 	print_int_vector(n, beta);
// 	Rprintf("Gamma:");
// 	print_int_vector(n, gamma);
// 	Rprintf("m: ");
// 	print_int_vector(1, m);
// 	Rprintf("k: ");
// 	print_int_vector(1, k);
// 
// 	Rprintf("LEFT MERGE!\n");
// 
// 	dictionary * tmp  = old_log_subperms;
// 	old_log_subperms = new_log_subperms;
// 	new_log_subperms = tmp;
// 
// 	sparse_reverse_lm(old_log_subperms, new_log_subperms, 4,  1,
// 				alpha, beta,  gamma,  log_factorials, n, m,
// 				k);
// 	//print_matrix(n+1, n+1, new_log_subperms);
// 	print_sparse_matrix(new_log_subperms);
// 	Rprintf("Vectors:\n Alpha: ");
// 	print_int_vector(n, alpha);
// 	Rprintf("Beta:");
// 	print_int_vector(n, beta);
// 	Rprintf("Gamma:");
// 	print_int_vector(n, gamma);
// 	Rprintf("m: ");
// 	print_int_vector(1, m);
// 	Rprintf("k: ");
// 	print_int_vector(1, k);
// 
// 
// 
// 
// 	tmp  = old_log_subperms;
// 	old_log_subperms = new_log_subperms;
// 	new_log_subperms = tmp;
// 
// 
// 	Rprintf("BOTTOM SPLIT!\n");
// 	Rprintf("old matrix:\n");
// 	print_sparse_matrix(old_log_subperms);
// 	Rprintf("new matrix before update:\n");
// 	print_sparse_matrix(new_log_subperms);
// 	sparse_reverse_bs(old_log_subperms, new_log_subperms, 1,  0,
// 				alpha, beta,  gamma,  log_factorials, n, m,
// 				k);
// 	Rprintf("RESULT:\n");
// 	print_sparse_matrix(new_log_subperms);
// 	Rprintf("Vectors:\n Alpha: ");
// 	print_int_vector(n, alpha);
// 	Rprintf("Beta:");
// 	print_int_vector(n, beta);
// 	Rprintf("Gamma:");
// 	print_int_vector(n, gamma);
// 	Rprintf("m: ");
// 	print_int_vector(1, m);
// 	Rprintf("k: ");
// 	print_int_vector(1, k);
// 
// 
// 	tmp  = old_log_subperms;
// 	old_log_subperms = new_log_subperms;
// 	new_log_subperms = tmp;
// 
// 
// 	Rprintf("TOP TRIM!\n");
// 	Rprintf("old matrix:\n");
// 	print_sparse_matrix(old_log_subperms);
// 	Rprintf("new matrix before update:\n");
// 	print_sparse_matrix( new_log_subperms);
// 	sparse_reverse_tt(old_log_subperms, new_log_subperms, 1,  0,
// 				alpha, beta,  gamma,  log_factorials, n, m,
// 				k);
// 	Rprintf("RESULT:\n");
// 	print_sparse_matrix( new_log_subperms);
// 	Rprintf("Vectors:\n Alpha: ");
// 	print_int_vector(n, alpha);
// 	Rprintf("Beta:");
// 	print_int_vector(n, beta);
// 	Rprintf("Gamma:");
// 	print_int_vector(n, gamma);
// 	Rprintf("m: ");
// 	print_int_vector(1, m);
// 	Rprintf("k: ");
// 	print_int_vector(1, k);
// 
// 
// 
// 
// 	tmp  = old_log_subperms;
// 	old_log_subperms = new_log_subperms;
// 	new_log_subperms = tmp;
// 
// 	
// 
// 	Rprintf("RIGHT MERGE!\n");
// 	Rprintf("old matrix:\n");
// 	print_sparse_matrix( old_log_subperms);
// 	Rprintf("new matrix before update:\n");
// 	print_sparse_matrix(new_log_subperms);
// 	sparse_reverse_rm(old_log_subperms, new_log_subperms, 1,  1,
// 				alpha, beta,  gamma,  log_factorials, n, m,
// 				k);
// 	Rprintf("RESULT:\n");
// 	print_sparse_matrix( new_log_subperms);
// 	Rprintf("Vectors:\n Alpha: ");
// 	print_int_vector(n, alpha);
// 	Rprintf("Beta:");
// 	print_int_vector(n, beta);
// 	Rprintf("Gamma:");
// 	print_int_vector(n, gamma);
// 	Rprintf("m: ");
// 	print_int_vector(1, m);
// 	Rprintf("k: ");
// 	print_int_vector(1, k);
// 
// 
// 
// 	tmp  = old_log_subperms;
// 	old_log_subperms = new_log_subperms;
// 	new_log_subperms = tmp;
// 
// 
// 
// 	Rprintf("TOP SPLIT!\n");
// 	//Rprintf("old matrix:\n");
// 	//print_sparse_matrix( old_log_subperms);
// 	//Rprintf("new matrix before update:\n");
// 	//print_matrix(n+1, n+1, new_log_subperms);
// 	sparse_reverse_ts(old_log_subperms, new_log_subperms, 1,  1,
// 				alpha, beta,  gamma,  log_factorials, n, m,
// 				k);
// 	Rprintf("RESULT:\n");
// 	print_sparse_matrix(new_log_subperms);
// 	Rprintf("Vectors:\n Alpha: ");
// 	print_int_vector(n, alpha);
// 	Rprintf("Beta:");
// 	print_int_vector(n, beta);
// 	Rprintf("Gamma:");
// 	print_int_vector(n, gamma);
// 	Rprintf("m: ");
// 	print_int_vector(1, m);
// 	Rprintf("k: ");
// 	print_int_vector(1, k);
// 
// 
// 
// 	tmp  = old_log_subperms;
// 	old_log_subperms = new_log_subperms;
// 	new_log_subperms = tmp;
// 
// 
// 
// 	Rprintf("LEFT MERGE!\n");
// 	//Rprintf("old matrix:\n");
// 	//print_matrix(n+1,n+1, old_log_subperms);
// 	//Rprintf("new matrix before update:\n");
// 	//print_matrix(n+1,n+1, new_log_subperms);
// 	sparse_reverse_lm(old_log_subperms, new_log_subperms, 1,  3,
// 				alpha, beta,  gamma,  log_factorials, n, m,
// 				k);
// 	Rprintf("RESULT:\n");
// 	print_sparse_matrix( new_log_subperms);
// 	Rprintf("Vectors:\n Alpha: ");
// 	print_int_vector(n, alpha);
// 	Rprintf("Beta:");
// 	print_int_vector(n, beta);
// 	Rprintf("Gamma:");
// 	print_int_vector(n, gamma);
// 	Rprintf("m: ");
// 	print_int_vector(1, m);
// 	Rprintf("k: ");
// 	print_int_vector(1, k);
// 
// 
// 	tmp  = old_log_subperms;
// 	old_log_subperms = new_log_subperms;
// 	new_log_subperms = tmp;
// 
// 	
// 
// 	Rprintf("BOTTOM SPLIT!\n");
// 	//Rprintf("old matrix:\n");
// 	//print_matrix(n+1, n+1, old_log_subperms);
// 	//Rprintf("new matrix before update:\n");
// 	//print_sparse_matrix( new_log_subperms);
// 	sparse_reverse_bs(old_log_subperms, new_log_subperms, 1,  3,
// 				alpha, beta,  gamma,  log_factorials, n, m,
// 				k);
// 	Rprintf("RESULT:\n");
// 	print_sparse_matrix(new_log_subperms);
// 	Rprintf("Vectors:\n Alpha: ");
// 	print_int_vector(n, alpha);
// 	Rprintf("Beta:");
// 	print_int_vector(n, beta);
// 	Rprintf("Gamma:");
// 	print_int_vector(n, gamma);
// 	Rprintf("m: ");
// 	print_int_vector(1, m);
// 	Rprintf("k: ");
// 	print_int_vector(1, k);
// 
// // 
// // 
// // 	// getting reduced logperms
// // 
// // 
// // 
// // 	alpha[0] = 5;
// // 	alpha[1] = 2;
// // 
// // 	beta[0] = 2;
// // 	gamma[0] =0;
// // 	*m = 2;
// // 	*k = 2;
// // 
// // 
// // 
// // 	Rprintf("REDUCED SUBPERMS\n");
// // 	get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
// // 				log_factorials, n,  m, k);
// // 	Rprintf("RESULT:\n");
// // 	print_matrix(n+1, n+1, new_log_subperms);
// // 
// // 
// // 
// 	// trying reduced!!
// 
// 	alpha[0] = 1;
// 	alpha[1] = 3;
// 	alpha[2] = 1;
// 	alpha[3] = 1;
// 	alpha[4] = 1;
// 
// 
// 	beta[0] = 0;
// 	beta[1] = 1;
// 	beta[2] = 2;
// 	beta[3] = 1;
// 
// 	gamma[0] =1;
// 	gamma[1] =1;
// 	gamma[2] =1;
// 	gamma[3] =0;
// 
// 	*m = 7;
// 	*k = 5;
// 
// 	int history_len = 0;
// 
// 	SEXP historySEXP = PROTECT(allocVector(INTSXP, 3*n));
// 	SEXP amount_historySEXP = PROTECT(allocVector(INTSXP, 6*n));
// 
// 	int * history = INTEGER(historySEXP);
// 	int * amount_history = INTEGER(amount_historySEXP);
// 
// 	memset(history, 0, sizeof(int)*3*n);
// 	memset(amount_history, 0, sizeof(int)*6*n);
// 
// 	Rprintf("REDUCING NOW\n");
// 	reduction(alpha,  beta,  gamma, m, n, k, history,
// 			   amount_history, &history_len);
// 
// 
// 
// 	Rprintf("REDUCED SUBPERMS\n");
// 	sparse_get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
// 				log_factorials, n,  m, k);
// 	Rprintf("RESULT:\n");
// 	print_sparse_matrix(new_log_subperms);
// 
// 	tmp  = old_log_subperms;
// 	old_log_subperms = new_log_subperms;
// 	new_log_subperms = tmp;
// 
// 
// 
// 	Rprintf("==========\nReverse reduction:\n==========\n");
// 	Rprintf("old = %d\n", old_log_subperms);
// 	Rprintf("new = %d\n", new_log_subperms);
// 	dictionary * the_log_subperms = sparse_reverse_reduction(old_log_subperms, new_log_subperms, alpha,
// 					   beta,  gamma, m,  n, k,  history,
// 			           amount_history, &history_len, log_factorials);
// 
// 
// 	for (int i = 0; i < history_len; ++i)
// 	{
// 		Rprintf("amounthistory[%d] = (%d,%d)\n", i, amount_history[2*i], amount_history[2*i+1]);
// 	}
// 
// 
// 	print_sparse_matrix(the_log_subperms);
// 	double logperm =  Csparse_log_sum_exp(the_log_subperms);
// 	Rprintf("logperm = %f\n", logperm);
// 
// 
// 	free_dictionary(new_log_subperms);
// 	free_dictionary(old_log_subperms);
// 	UNPROTECT(9);
// 	
// 
// 
// 
// 	return(mSEXP);
// 
// }