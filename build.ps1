# Build the code to produce an executable
g++ -Wall main.cpp mat.cpp -o out
# Run the executable
./out
#
Remove-Item -Force -ErrorAction SilentlyContinue out.exe
