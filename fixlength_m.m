function fixlength2(s1,s2,len,indent,fid)
%FIXLENGTH Returns a string which has been divided up into < LEN
%character chuncks with the '...' line divide appended at the end
%of each chunck.
%   FIXLENGTH(S1,S2,L) is string S1 with string S2 used as break points into
%   chuncks less than length L.

%Eric Westervelt
%5/31/00
%1/11/01 - updated to use more than one dividing string
%4/29/01 - updated to allow for an indent
%3/24/05 - bjm - non-iterative operation, much faster for large files

	A=[];
	for c = 1:length(s2)
		A = [A findstr(s1,s2(c))];
		A = sort(A);
	end
	
	if (length(A)>0)
		
		dA = [diff(A)];
		B = A;
		
		runlen = A(1);
		
		for i=1:1:length(dA)
			if (runlen + dA(i) <= len)
				runlen = runlen + dA(i);
				B(i) = 0;
			else
				runlen = 0;
			end
		end
		B(end) = 0; % don't make a newline after the last entry
		
		B=[0 B(B>0) length(s1)];
		
		for i=1:1:length(B)-1
			if (i==1)	
				fprintf(fid,'%s',s1(B(i)+1:B(i+1)));
			else
				fprintf(fid,'...\n%s%s',indent,s1(B(i)+1:B(i+1)));
			end
		end
		
	else
		fprintf(fid,'%s',s1);
	end