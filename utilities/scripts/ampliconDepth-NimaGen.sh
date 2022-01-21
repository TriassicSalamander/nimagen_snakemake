#Developed by: Dr. Sreenu Vattipally at University of Glasgow

samtools view -h $2 |\
	awk -F '\t' '{
		if(NR==FNR) {
			pStart[NR]=$3; pEnd[NR]=$6; pCount++;
		}
		else {
			if ($1!~/^@/&&!and($2,0x4)&&!and($2,0x800)) {
				end=0;
				start=$4; gsub(/[A-Z]/," & ",$6); n=split($6,a," ");
				for(i=0;i<=n;i++)
					if(a[i]=="M")
						end+=a[i-1];
				for(i=1;i<=pCount;i++) {
					if(start>=pStart[i]&&((start+end)<=pEnd[i]) && end > 30) {
						Primer[i]++;
					}
				}
			}
		}
	};
	END{
		for(i=1;i<=pCount;i++)
			print i"\t"pStart[i]"\t"pEnd[i]"\t"Primer[i]+0;
		}' $1 -
