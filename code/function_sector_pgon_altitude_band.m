
function sectors_pgon = function_sector_pgon_altitude_band(sector_data, altitude_band)

% Altitude band: FL

N = length(sector_data); %Number of sectors

sectors_pgon = cell(size(sector_data));

for i = 1:N
   sector_i = sector_data{i}; 
   
   % Select those in the altitude band (+-500ft)
   aux = [sector_i.properties]; aux2 = [aux.UPPER_LIMIT_VALUE]; upper_FL = aux2';
   aux = [sector_i.properties]; aux2 = [aux.LOWER_LIMIT_VALUE]; lower_FL = aux2';
   A = altitude_band*ones(size(lower_FL)); sec_index = (lower_FL<A)&(upper_FL>A); 
   sector_i(~sec_index)=[];
   
   Ni = length(sector_i); % Number of subsectors in sector i in the altitude band
   
   if Ni > 0
       lon1 = sector_i(1).geometry.coordinates(:,:,1);
       lat1 = sector_i(1).geometry.coordinates(:,:,2);
       
       pgon = polyshape(lon1, lat1);
       if ~ispolycw(lon1, lat1)
        [lon1, lat1] = poly2cw(lon1, lat1);
       end
       
       if Ni>1
           for j = 2:Ni
               lon2 = sector_i(j).geometry.coordinates(:,:,1);
               lat2 = sector_i(j).geometry.coordinates(:,:,2);
               
               if ~ispolycw(lon2, lat2)
                   [lon2, lat2] = poly2cw(lon2, lat2);
               end
                   [lon1, lat1] = polybool('union', lon1, lat1, lon2, lat2);
                   pgon = polyshape(lon1, lat1);               
           end
       end
       
       % Save polygon
       sectors_pgon{i} = pgon;
   end
end

end