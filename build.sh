dir="${PREFIX}/share"
mkdir -p "${dir}"
mkdir -p "${PREFIX}/bin"
cp recognizer.py "${dir}"
cp -r resources/* "${dir}"
chmod +x "${dir}/recognizer.py"
ln -s "${dir}/recognizer.py" "${PREFIX}/bin/recognizer"