g++ -std=c++17 main.cpp -o monotone
if [ $? -eq 0 ]; then
    ./monotone
    python3 plot.py
else
    echo "Compilation failed!"
fi