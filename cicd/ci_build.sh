PREFIX="/opt/conda"
dir="${PREFIX}/share"
mkdir -p "${dir}" "${PREFIX}/bin"
cp reCOGnizer/recognizer.py reCOGnizer/resources/* "${dir}"
chmod +x "${dir}/recognizer.py"
ln -s "${dir}/recognizer.py" "${PREFIX}/bin/recognizer"