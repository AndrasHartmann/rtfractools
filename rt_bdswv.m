function result=rt_bdswv(data,p)
%function result=rt_bdswv(data,p)
%Real-time bridge detrended scaled window variance analysis with window size 6^p
%
%Not really real-time though, because we feed in the entire input, but we will simulate as if we would be real-time
%
%    RT_Fractools: Real-time fractal signal processing in the time domain
%    Copyright (C) {2013} Andras Hartmann <andras.hartmann@gmail.com>
%
%    Please, if you use this program for any scientific purpose, cite 
%    Hartmann, A., Mukli, P., Nagy, Z., Kocsis, L., Hermán, P., & Eke, A. (2013).
%    Real-time fractal signal processing in the time domain. Physica A:
%    Statistical Mechanics and Its Applications, 392(1), 89–102.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%helper functions just for modularity (keep it simple)
%this one is for accumulating x and x^2 into sx and sx2, starting a new window
function x1 = bdswv_helper_1()
    xstart(j)=x;
    if (j==1)
	sx(j)=x;
	sx2(j)=x.*x;
    end;
    sxi(j)=x.*nn;
end;

%this one calculates the result at the end of the sliding window lgsd is the log of SD
function x2 = bdswv_helper_res()
    for k=1:n;
	lgsd(k)=log(sum(ssxd(k,:)./2^(p-k-used1+1)));
	lgn(k)=log(2^(k+used1-1));
    end;
    pp=polyfit(lgn,lgsd,1);
    result(value,2) = i;
    result(value,1)=pp(1) ./ 2;
    value = value + 1;
end;

%otherwise: this is the normal functionality accumulating within the window
function x3 = bdswv_helper_otherwise()
    if (j==1)
	sx(j)=sx(j)+x;
	sx2(j)=sx2(j)+x.*x;
    end;
    sxi(j)=sxi(j)+x.*nn;
end;

%this is to accumulate the result at the end of any window, not neccessary end of the sliding window
function x4 = bdswv_helper_0()
    mindex(j) = mindex(j) + 1;
    if (j==1)
	sx(j)=sx(j)+x;
	sx2(j)=sx2(j)+x.*x;
	ssx(j,mindex(j)) = sx(j);
	ssx2(j,mindex(j)) = sx2(j);
    end;
    if (mod(mindex(j),2)==0)
	ssx(j+1, mindex(j)./2) = ssx(j,mindex(j)) + ssx(j,mindex(j)-1);
	ssx2(j+1, mindex(j)./2) = ssx2(j,mindex(j)) + ssx2(j,mindex(j)-1);
    end;
    sxi(j)=sxi(j)+x.*nn;
    m=(x-xstart(j))./(nn-1); %calculating the bridge 

    %And here comes the magic, the nasty derived equation
    ssxd(j,mindex(j))=(ssx2(j,mindex(j))-2.*m.*sxi(j)+m.*(nn+1).*ssx(j,mindex(j))-ssx(j,mindex(j)).^2./nn+(m.^2.*nn.*(nn.^2-1))./12)./(nn-1);

if ssxd(j,mindex(j))<0
	ssxd(j,mindex(j)) = 0;
	%warning('Results might be compromised because of numerical instability');
end

end;

%N=data length
N=2^p;
% determining the used range (from used1 to used2)
ignored = [7 2 1; 8 2 1; 9 2 1; 10 2 2; 11 2 2; 12 2 4; 13 3 5; 14 3 5; 15 3 4; 16 3 6; 17 3 6];
% which bins to ignore: length, smallest bin size, largest bin size;
i=1;
while (p>ignored(i,1))
	i=i+1;
end;
smallest=ignored(i,2);
largest=ignored(i,3);
% smallest: ignored these smallest bin size;
% largest: ignored these largest bin size;
used1=smallest+1;
used2=p-largest;
%n = number of used points
n=used2-used1+1;

%summed values
sx=zeros(n);
sx2=zeros(n);
sxi=zeros(n);
ssx=zeros(n,n);
ssx2=zeros(n,n);
ssxd=zeros(n,n);
xstart=zeros(n);
mindex=zeros(n);
value = 1;

%in case you would like to measure the time (don't forget to toc too)
%tic

x = 0;
mindex=zeros(n);
for i=1:length(data);
    x=x+data(i);
    for j=1:n
	c = mod(i,2.^(j+used1-1)); %period
	nn = mod(i-1,2.^(j+used1-1))+1; %n
        switch c;
            case 1
		bdswv_helper_1;

            case 0
		bdswv_helper_0;
		if (i>=N && j == n)
		    bdswv_helper_res;
		end;

	   otherwise
		bdswv_helper_otherwise;

	end;
    end;
    if (mod(i,N)==0)
	mindex=zeros(n);
    end;

end;

%in case you would like to measure the time
%toc
%etime(clock,t0)
end

