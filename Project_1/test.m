fid = fopen('jj.txt','w');
ao=length(a);
for i=1:ao
    for j=1:ao
        if j==ao
            fprintf(fid,'%g \\\\',a(i,j));
        else
            fprintf(fid,'%g &',a(i,j));
        end
    end
end