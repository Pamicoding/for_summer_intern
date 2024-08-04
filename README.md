## h3dd steps
1. gamma2h3dd.py     
把gamma_picks.csv, gamma_events.csv轉成h3dd的格式    
2. station_append.py    
在測站資訊後面補上一些常數    
3. tomops_H14    
用來重定位的3維速度模型 (Hsu et al. 2020)     
4. h3dd.inp    
h3dd的參數檔，把前三項依序填上即可    
5. ./h3dd <h3dd.inp     
h3dd的主程式，直接在終端機執行即可，執行完就會產生下面兩個檔案         
如果不能用的話，我們需要自己compile fortran的code，進入到./H3DD/src，在終端機打make去compile程式     
```
make
```
* .dout    
重定位完的結果    
* .hout    
event_list    
