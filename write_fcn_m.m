function write_fcn(fcn_name,arguments,replace_list,list)

disp(['writing ',fcn_name]); 	% open file and write first line
fid = fopen(fcn_name,'w');

fprintf(fid,'function [');
for item = 1:1:size(list,1)
	currentname = list{item,2};
	if (item > 1) fprintf(fid,','); end
	fprintf(fid,'%s',currentname);
end

[path,name,ext] = fileparts(fcn_name);
fprintf(fid,'] = %s(',name);

for item = 1:1:size(arguments,2) % write arguments
	currentname = arguments{1,item};
	if (item > 1) fprintf(fid,','); end
	fprintf(fid,'%s',currentname);
end
fprintf(fid,')\n\n');

fprintf(fid,'  model_params;\n\n');

for item = 1:1:size(list,1) % write variables to file
	currentvar = list{item,1};
	currentname = list{item,2};
	[n,m]=size(currentvar);
	for i=1:n
		for j=1:m
			Temp0=currentvar(i,j);
			Temp1=replace(char(Temp0),replace_list);
			Temp2=['  ',currentname,'(',num2str(i),',',num2str(j),')=',Temp1,';'];
			fixlength_m(Temp2,'*+',100,'         ',fid);
			fprintf(fid,'\n');
		end
	end
	fprintf(fid,'\n%s',' ');
end

status = fclose(fid);
disp(' - done');
%=============================================

function str = replace(str,replace_list)
	
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

	for i=size(replace_list,1):-1:1
		orig_str = char(replace_list(i,1));
		new_str = char(replace_list(i,2));
		str = strrep(str,orig_str,new_str);
	end
