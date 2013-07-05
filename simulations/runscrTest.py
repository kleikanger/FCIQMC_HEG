import runcode
a  = runcode.runcode()

a.i_numclosedshells = 2
a.dS = 5
a.i_num_loops = 100000
a.i_initiatorlimit = 3
a.d_rs = 1

a.i_numshells_basis = 3
for i in range (2,10,4):
	a.i_limit_nw = i*1000
	a.s_ofilename = '10pt_2D_rs1_4shell_i3_nw%sk'%i
	a.runcode()

a.i_numshells_basis = 4
for i in range (2,14,4):
	a.i_limit_nw = i*1000
	a.s_ofilename = '10pt_2D_rs1_5shell_i4_nw%sk'%i
	a.runcode()

a.i_numshells_basis = 5
for i in range (2,18,4):
	a.i_limit_nw = i*1000
	a.s_ofilename = '10pt_rs1_2D_6shell_i5_nw%sk'%i
	a.runcode()

a.i_numshells_basis = 6
for i in range (2,18,4):
	a.i_limit_nw = i*1000
	a.s_ofilename = '10pt_rs1_2D_6shell_i6_nw%sk'%i
	a.runcode()
