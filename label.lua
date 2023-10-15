-----------------------------------------------------------------------------------------------------------------------------------
--
-- Author: Rodrigo A. Moreira (C) 2023
-- https://orcid.org/0000-0002-7605-8722
-- LICENSE: CC BY-NC-ND 4.0 (https://creativecommons.org/licenses/by-nc-nd/4.0/)
--
----------------------------------------------------------------------------------------------------------------------------------

local letters = {'G','H','I','B','E','T','C'}
for k=1,#letters do letters[letters[k]] = k end

--[[ FROM STRIDE DOCUMENTATION

**) One-letter secondary structure code is nearly the same as used  in
           DSSP [2] (see Frishman and Argos [1] for details):

           H        Alpha helix
           G        3-10 helix
           I        PI-helix
           E        Extended conformation
           B or b   Isolated bridge
           T        Turn <-- hydrogen bond turn
-----------------------------------------------------
           C        STRIDE Coil (or none of the above)
           S        DSSP bend
--]]

local PDB = assert(arg[1])
local NETWORK = 'network_backboneRE_heavy_gt2'

local map = function(w)
	l3 = {'ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','SEC','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP','!!!'}
	l1= {'R','H','K','D','E','S','T','N','Q','C','U','G','P','A','V','I','L','M','F','Y','W','!'}
	assert(#l3 == #l1)
	for k=1,#l3 do
		if l3[k] == w or l1[k] == w then return l3[k],l1[k] end
	end
end

local resdb = function()
	local data = {}
	data.str = {}
	data.res = {}
	local count = 0
	return function(res,str)
		if not res then
			--print('resdb rescount',count,res,str)
			return data
		end
		assert(#res == 3)
		assert(#str == 1)
		if not letters[str] then str = 'C' end
		count = count + 1
		data.res[count] = res
		data.str[count] = str
	end
end

local dssp = false
do
-- DSSP example of output
--[[
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
    1    1 A M              0   0  158      0, 0.0   202,-0.2     0, 0.0     2,-0.2   0.000 360.0 360.0 360.0-134.4    9.7   30.8   70.3
    2    2 A E  P>    -     0   0   77    200,-2.0   200,-0.4     1,-0.1     2,-0.1  -0.518 360.0-113.7 -83.5 153.4    7.9   34.1   69.9
  444  444 A H    <         0   0    6     -3,-2.8    -4,-0.0     3,-0.1   159,-0.0  -0.700 360.0 360.0 -84.8 140.2   -5.0   55.0  112.8
  445  445 A E              0   0  198     -2,-0.3     3,-0.1     0, 0.0    -3,-0.0  -0.571 360.0 360.0  71.7 360.0   -3.7   52.8  115.6
  446        !              0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0
  447  469 A Q              0   0  185      2,-0.0    -3,-0.1     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 141.2    3.1   54.3  114.9
  448  470 A I        -     0   0   12     -3,-0.1    -4,-0.0     1,-0.1     0, 0.0  -0.203 360.0-126.3 -49.7 122.3    1.5   57.6  116.0
--]]
	local res = resdb()
	local readline = false
	local count = 0
	for line in io.lines(string.format("%s.dssp",PDB)) do
		if readline then --> is a residue
			local residue = line:sub(14,14)
			local terminal = line:sub(15,15)
			--if residue == '!' then print('WARNING:',count,'->',line) end
			if residue ~= '!' then
				count = count+1
				residue = residue:upper() ~= residue and 'CYS' or residue
				residue = map(residue)
				local letter = string.upper(line:sub(17,17))
				res(residue,letter)
			end
		end
		if line:find("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC") then readline = true end
	end
	dssp = res()
end

local stride = false
do
-- STRIDE example output
--[[
REM                                                                        1MJH
REM  --------------- Detailed secondary structure assignment-------------  1MJH
REM                                                                        1MJH
REM  |---Residue---|    |--Structure--|   |-Phi-|   |-Psi-|  |-Area-|      1MJH
ASG  VAL A    3    1    C          Coil    360.00    131.52      82.9      1MJH
ASG  MET A    4    2    C          Coil   -136.42    126.23     118.7      1MJH
ASG  TYR A    5    3    T          Turn    -89.26     75.08      25.9      1MJH
ASG  LYS A    6    4    T          Turn    -87.37    -31.13     112.8      1MJH
ASG  LYS A    7    5    E        Strand   -130.53    121.59      49.7      1MJH
ASG  ILE A    8    6    E        Strand   -111.05    138.26       0.6      1MJH
ASG  LEU A    9    7    E        Strand   -101.13    120.63       0.2      1MJH
ASG  TYR A   10    8    E        Strand   -115.42     94.78       2.8      1MJH
ASG  PRO A   11    9    E        Strand    -85.43    122.97      12.1      1MJH
ASG  THR A   12   10    C          Coil   -121.22    141.81       8.1      1MJH
ASG  ASP A   13   11    C          Coil   -108.11     11.06      54.2      1MJH
ASG  PHE A   14   12    C          Coil     63.69     23.94      34.1      1MJH
ASG  SER A   15   13    C          Coil    -85.04    170.71       2.4      1MJH
ASG  GLU A   16   14    H    AlphaHelix    -68.26    -33.20     111.1      1MJH
ASG  THR A   17   15    H    AlphaHelix    -74.71    -34.01       3.0      1MJH
ASG  ALA A   18   16    H    AlphaHelix    -60.70    -32.06       2.8      1MJH
ASG  GLU A   19   17    H    AlphaHelix    -68.76    -32.66     112.7      1MJH
ASG  ILE A   20   18    H    AlphaHelix    -64.80    -43.93      38.9      1MJH
ASG  ALA A   21   19    H    AlphaHelix    -62.36    -30.38       0.0      1MJH
ASG  LEU A   22   20    H    AlphaHelix    -60.69    -31.50      12.0      1MJH
ASG  LYS A   23   21    H    AlphaHelix    -64.55    -36.23     125.7      1MJH
--]]
	local res = resdb()
	for line in io.lines(string.format('%s.stride',PDB)) do if line:sub(1,3) == 'ASG' then
		local residue = map(line:sub(6,8))
		local letter = string.upper(line:sub(25,25))
		res(residue,letter)
	end end
	stride = res()
end

-- Check residues sequence
for k=1,math.min(#(dssp.res),#(stride.res)) do
	if dssp.res[k] ~= '!!!' and dssp.res[k] ~= stride.res[k] then error(string.format('%s %s DSSP.RES[%d] = %s BUT STRIDE.RES[%d] = %s',PDB,NETWORK,k,dssp.res[k],k,stride.res[k])) end
end
local rescount = #(dssp.res)

local split = function(str)
	local t = {}
	for w in str:gmatch('%S+') do t[#t+1] = w end
	return t
end

local oldcut = 0.0
local words = false
for line in io.lines(string.format("%s.%s.residues_curvature",PDB,NETWORK)) do
	local words = split(line)
	local cut = words[1]

	oldwords = oldwords or words
	local oldcut = oldwords[1]

	do
	local c2 = math.floor(tonumber(cut)*10)
	local c1 = math.floor(tonumber(oldcut)*10)+1
	while c1 < c2 do
		for k=2,rescount+1 do
			io.write(string.format('%.3f %s %d %s %s %s\n',c1/10,dssp.res[k-1],k-1,dssp.str[k-1],stride.str[k-1],oldwords[k]))
		end
		c1 = c1 + 1
	end
	end
	for k=2,rescount+1 do
		io.write(string.format('%s %s %d %s %s %s\n',cut,dssp.res[k-1],k-1,dssp.str[k-1],stride.str[k-1],words[k]))
	end
	
	oldwords = words
end


