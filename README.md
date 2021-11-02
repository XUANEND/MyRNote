# MyRNote

有时候，推送时会出现下面的问题：

failed to connect to github.com port 443

解决办法就是设置代理

git 代理设置（本机）：

` git config --global http.proxy http://127.0.0.1:1081`   
` git config --global https.proxy http://127.0.0.1:1081`

取消代理：

` git config --global --unset http.proxy`  
` git config --global --unset https.proxy`

