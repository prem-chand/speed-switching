function write_fcn_c(fcn_name,inputs,replace_list,outputs)

% Writes a function in C
%
%   write_fcn_c(fcn_name, inputs, replace_list, outputs)
%
%   fcn_name    A character string denoting the name of the function we are
%       writing
%   inputs      A two column list of the input parameters with the fist
%       column being the data type and the second the variable name.
%   replace_list    A two column matrix of strings with the first column
%       containing the string to replace with the string in the second
%       column.
%   outputs     A three column list of the output parameters with the first
%       column containing the variable, the second its c variable type, and
%       the third its name.

%--------------------------------------------------------------------------
% process the longest strings first to keep from accidentally replacing substrings of longer strings
% perform a bubble sort on 'replace_list' according to decending length of the first element

for loop=1:1:size(replace_list,1)
	for i=1:1:size(replace_list,1)-1
		if (length(replace_list{i,1}) < length(replace_list{i+1,1}))
			temp = replace_list(i,:);
			replace_list(i,:) = replace_list(i+1,:);
			replace_list(i+1,:) = temp;
		end
	end
end

disp(['writing ',fcn_name]); 	% open file and write first line
fid = fopen(fcn_name,'w');

fprintf(fid,'void ');

[path,name,ext,ver] = fileparts(fcn_name);
fprintf(fid,'%s(',name);

for item = 1:1:length(inputs) % write inputs
	currentname = inputs{item};
	if (item > 1) fprintf(fid,', '); end	
	fprintf(fid,'%s %s','dMatrix', currentname);
end

for item = 1:1:size(outputs,1)  % write outputs
	currentname = outputs{item,2};
	fprintf(fid,', ');
	fprintf(fid,'%s *%s','dMatrix', currentname);
end

%fprintf(fid,',%s %s','struct PARAMS', '*params');

fprintf(fid,')\n{\n\n');

for i=1:1:length(replace_list)
	fprintf(fid,'  double %s = %s;\n',char(replace_list(i,1)),char(replace_list(i,2)));
end
fprintf(fid,'\n');

fid_tem = fopen('temp.tmp','w');
declare_list = {};
s='';
for item = 1:1:size(outputs,1) % write variables to file
	currentvar = outputs{item,1};
	currentname = outputs{item,2};
	[n,m]=size(currentvar);

	fprintf(fid,'  dMatrix %s_tem(%i,%i);\n',currentname,n,m);
	fprintf(fid,' *%s = %s_tem;\n\n',currentname,currentname);

	list_tem = cc(fid_tem,currentvar,currentname);
	declare_list = [declare_list;list_tem];
end

fclose(fid_tem);

%=======================================
% remove duplicates in declaration list
%=======================================

declare_list=sort(declare_list);
unique=ones(size(declare_list,1),1);
for i=1:1:size(declare_list,1)-1
	if (strcmp(char(declare_list(i)),char(declare_list(i+1))))
		unique(i) = 0;
	end
end
declare_list = declare_list(unique == 1);

for i=1:1:size(declare_list,1)
	fprintf(fid,'  double %s;\n',char(declare_list(i)));
end

fprintf(fid,'\n');

%====================
% print data to file
%====================

fid_tem = fopen('temp.tmp','r');

while(1)
	line = fgetl(fid_tem);
	if ~ischar(line), break, end
	fprintf(fid,'%s\n',line);
end

fprintf(fid, '  return;\n}\n\n');
fclose(fid_tem);
fclose(fid);
delete('temp.tmp');
disp(' - done');

%=============================================

function [str] = replace(str,replace_list) 
	
  replaced = zeros(size(str));
	for i=size(replace_list,1):-1:1
		orig_str = char(replace_list(i,1));
		new_str = char(replace_list(i,2));

		t = findstr(str,orig_str);

		while(any(t))
			ok = 1;

			for i=t(1):1:t(1)+length(orig_str)-1
				if (replaced(i) == 1)
					ok = 0;
				end
			end

			if (ok == 1)
				str = [str(1:t(1)-1), new_str, str(t(1)+length(orig_str):end)];
				replaced = [replaced(1:t(1)-1), ones(size(new_str)),replaced(t(1)+length(orig_str):end)];			
				t = t(2:end) + ones(size(t(2:end)))*(length(new_str)-length(orig_str));
			else
				t = t(2:end);
			end
		end
	end

%==================================================================

function [declare_list] = cc(fid,var,name)

str = char(var);
if (length(findstr(str,'matrix'))==0)
	str = ['matrix([[', str, ']])'];
end

c = optimized_ccode(name,str);

% string replacements necessary for correct syntax
c = strrep(c,' ','');
c = strrep(c,'~','');
c = strrep(c,'][',',');
c = strrep(c,'[','(');
c = strrep(c,']',')');
c = strrep(c,name,['(*',name,')']);

I = findstr(c,';');
I = [0 I];

declare_list = {};
for i=1:1:length(I)-1
	if (length(findstr(c(I(i)+1:I(i+1)),name))==0)
		fprintf(fid,'  %s\n',c(I(i)+1:I(i+1)));
		tem=findstr(c(I(i)+1:I(i+1)),'=');
		declare_list=[declare_list;{c(I(i)+1:I(i)+tem-1)}];
	else
		fprintf(fid,'  %s\n',c(I(i)+1:I(i+1)));
	end
end

fprintf(fid,'\n\n');