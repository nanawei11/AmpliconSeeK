ss=$1
find . -name *_ask_stats.csv |grep $ss |sort |xargs grep circular |grep -v ": 0$"
