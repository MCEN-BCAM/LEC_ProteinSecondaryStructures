-----------------------------------------------------------------------------------------------------------------------------------
--
-- Author: Rodrigo A. Moreira (C) 2023
-- https://orcid.org/0000-0002-7605-8722
-- LICENSE: CC BY-NC-ND 4.0 (https://creativecommons.org/licenses/by-nc-nd/4.0/)
--
----------------------------------------------------------------------------------------------------------------------------------

local PDB = arg[1]

local pdbatom = function(line)
--[[
ATOM     10 HG21 THR     1      55.970  81.180  25.060  1.00  0.00
ATOM      1  C1  BGLCM   1     -11.173   1.135-257.897  1.00  1.00      M0   C
--]]
	local mt = {
	__index = function(t,k)
		local line = t[1]
		local s = function(i,j) return string.upper(line:sub(i,j):gsub('%s+','')) end
		local i = function(i,j) return math.tointeger(line:sub(i,j)) end
		local f = function(i,j) return tonumber(line:sub(i,j)) end
		if k == 1  or k == 'serial'   then return i(7,11)  end
		if k == 2  or k == 'atomname' then return s(13,16) end
		if k == 3  or k == 'altloc'   then return s(17,17) end
		if k == 4  or k == 'resname'  then return s(18,21) end
		if k == 5  or k == 'chain'    then return s(22,22) end
		if k == 6  or k == 'resseq'   then return i(23,26) end
		if k == 7  or k == 'code'     then return f(27,27) end
		if k == 8  or k == 'x'        then return f(31,38) end
		if k == 9  or k == 'y'        then return f(39,46) end
		if k == 10 or k == 'z'        then return f(47,54) end
		if k == 11 or k == 'occ'      then return f(55,60) end
		if k == 12 or k == 'temp'     then return f(61,66) end
		if k == 13 or k == 'seg'      then return s(73,76) end
		if k == 14 or k == 'symb'     then return s(77,78) end
		if k == 15 or k == 'charge'   then return s(79,80) end
	end,
	__newindex = function(t,k,v) print('WARNING: BETTER NOT TO CHANGE ATOM PROPERTIES!') end,
	__len = function(t)
		local count = 0
		for k=1,15 do if t[k] then count = count + 1 end end
		return count
	end
	}
	assert(type(line) == 'string' and line:sub(1,4) == 'ATOM')
	local ret = setmetatable({line},mt)
	assert(#ret >= 10)
	return ret
end

local nodes = {}

local flag = false
for line in io.lines(PDB) do
	--if line:sub(1,4) ~= 'ATOM' and flag then break end -- select only first constinuous segment of each pdb
	if line:sub(1,4) == 'ATOM' then
		flag = true
		local t = pdbatom(line)
		local resname = t.resname
		if (resname == "CL" or resname == 'NA' or resname == 'SOL') then goto fim end
		local atomname = t.atomname
		if atomname:sub(1,1) == 'H' then  goto fim end
		local n = {'N','H','CA','HA','C'}
		for k=1,#n do if n[k] == atomname then goto add end end	
		goto fim
		::add::
		do
			local resid = string.format('%s%s%s',t.resname,t.chain=='' and ' ' or t.chain,t.resseq)
			nodes[#nodes + 1] = {t.x,t.y,t.z,resid,atomname,t.chain=='' and ' ' or t.chain,t.resseq}
		end
		::fim::
	end
end

local dist = function(a,b)
	local x1,y1,z1 = a[1],a[2],a[3]
	local x2,y2,z2 = b[1],b[2],b[3]
	local dx,dy,dz = x1-x2,y1-y2,z1-z2
	return math.sqrt(dx*dx + dy*dy + dz*dz)
end

data = {}
buffer = {}

local nnodes = #nodes
for k=1,#nodes do
for v=k+1,#nodes do
	local d = dist(nodes[k],nodes[v])
	if d>2.0 and d<7.8 then
		local resid1  = nodes[k][4]
		local resid2  = nodes[v][4]
		local atomname1 = nodes[k][5]
		local atomname2 = nodes[v][5]
		local chain1 = nodes[k][6]
		local chain2 = nodes[v][6]
		local resseq1 = nodes[k][7]
		local resseq2 = nodes[v][7]
		data[#data+1] = {k,v,d}
	end
end
end

io.write(nnodes,'\n')
for i=1,#data do
	local k,v,d = data[i][1],data[i][2],data[i][3]
	io.write(string.format("%d %d %f\n",k,v,d))
end











