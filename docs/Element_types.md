## Finite element types

|Name|IP location|
|:-:|:-:|
|Truss2n1ip|![This is an image](./assets/images/2n1ip.svg)</br>ip1 (r = 0)|
|Beam2n1ip|![This is an image](./assets/images/2n1ip.svg)</br>ip1 (r = 0)|
|Mem4n4ip|![This is an image](./assets/images/4n4ip.svg)</br>ip1 (r = 0.577, s = 0.577)</br>ip2 (r = -0.577, s = 0.577)</br>ip3 (r = -0.577, s = -0.577) </br>ip4 (r = 0.577, s = -0.577)|
|Plate4n4ipT|![This is an image](./assets/images/4n4ip.svg)</br>For tension and bending strain:</br>ip1 (r = 0.577, s = 0.577)</br>ip2 (r = -0.577, s = 0.577)</br>ip3 (r = -0.577, s = -0.577)</br>ip4 (r = 0.577, s = -0.577)</br>![This is an image](./assets/images/4n4ip_shear.svg)</br>For shear strain:</br>ipA (x = x1 - x2 - x3 + x4)</br>ipA (y = y1 - y2 - y3 + y4)</br>ipB (x = x1 - x2 + x3 - x4)</br>ipB (y = y1 - y2 + y3 - y4)</br>ipC (x = x1 + x2 - x3 - x4)</br>ipC (y = y1 + y2 - y3 - y4)|
