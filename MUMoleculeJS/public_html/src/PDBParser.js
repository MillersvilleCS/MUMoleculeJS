(function(window) {
    'use strict';

    var PDBParser = function( ) {

    };

    PDBParser.parse = function(text) {

        var atoms = [];
        var bonds = [];
        var bondHash = []; // used to check for multiple bonds
        var atomSerialMap = new Map();
        var protein = {sheet: [], helix: [], biomtChains: '', biomtMatrices: [], symMat: [], pdbID: '', title: ''};
        ;

        var lines = text.split('\n');

        for (var i = 0; i < lines.length; ++i) {

            var line = lines[i];

            if (line.substr(0, 4) === 'ATOM' || line.substr(0, 6) === 'HETATM') {

                parseAtom(line, atoms, atomSerialMap);

            } else if (line.substr(0, 6) === 'CONECT') {

                parseConnect(line, bonds, bondHash, atomSerialMap);

            } else if (line.substr(0, 5) === 'HELIX') {

                parseHelix(line, protein);

            } else if (line.substr(0, 6) === 'CRYST1') {

                protein.a = parseFloat(line.substr(6, 9));
                protein.b = parseFloat(line.substr(15, 9));
                protein.c = parseFloat(line.substr(24, 9));
                protein.alpha = parseFloat(line.substr(33, 7));
                protein.beta = parseFloat(line.substr(40, 7));
                protein.gamma = parseFloat(line.substr(47, 7));
                protein.spacegroup = line.substr(55, 11);
                defineCell(protein);

            } else if (line.substr(0, 6) === 'REMARK') {

                var type = parseInt(line.substr(7, 3));

                if (type === 290 && line.substr(13, 5) === 'SMTRY') {

                    var n = parseInt(line[18]) - 1;
                    var m = parseInt(line.substr(21, 2));

                    if (!protein.symMat[m])
                        protein.symMat[m] = new THREE.Matrix4().identity();

                    protein.symMat[m].elements[n] = parseFloat(line.substr(24, 9));
                    protein.symMat[m].elements[n + 4] = parseFloat(line.substr(34, 9));
                    protein.symMat[m].elements[n + 8] = parseFloat(line.substr(44, 9));
                    protein.symMat[m].elements[n + 12] = parseFloat(line.substr(54, 10));

                } else if (type === 350 && line.substr(13, 5) === 'BIOMT') {

                    var n = parseInt(line[18]) - 1;
                    var m = parseInt(line.substr(21, 2));

                    if (!protein.biomtMatrices[m])
                        protein.biomtMatrices[m] = new THREE.Matrix4().identity();

                    protein.biomtMatrices[m].elements[n] = parseFloat(line.substr(24, 9));
                    protein.biomtMatrices[m].elements[n + 4] = parseFloat(line.substr(34, 9));
                    protein.biomtMatrices[m].elements[n + 8] = parseFloat(line.substr(44, 9));
                    protein.biomtMatrices[m].elements[n + 12] = parseFloat(line.substr(54, 10));

                } else if (type === 350 && line.substr(11, 11) === 'BIOMOLECULE') {

                    protein.biomtMatrices = [];
                    protein.biomtChains = '';

                } else if (type === 350 && line.substr(34, 6) === 'CHAINS') {

                    protein.biomtChains += line.substr(41, 40);

                }
            } else if (line.substr(0, 6) === 'HEADER') {

                protein.pdbID = line.substr(62, 4);

            } else if (line.substr(0, 5) === 'TITLE ') {

                if (!protein.title)
                    protein.title = "";

                protein.title += line.substr(10, 70) + "\n"; // CHECK: why 60 is not enough???

            }
        }

        for (var i = 0; i < atoms.length; i++) {
            
            var atom = atoms[i];

            var found = false;
            // MEMO: Can start chain and end chain differ?
            for (var j = 0; j < protein.sheet.length; j++) {
                if (atom.chain != protein.sheet[j][0])
                    continue;
                if (atom.resi < protein.sheet[j][1])
                    continue;
                if (atom.resi > protein.sheet[j][3])
                    continue;
                atom.ss = 's';
                if (atom.resi == protein.sheet[j][1])
                    atom.ssbegin = true;
                if (atom.resi == protein.sheet[j][3])
                    atom.ssend = true;
            }
            for (j = 0; j < protein.helix.length; j++) {
                if (atom.chain != protein.helix[j][0])
                    continue;
                if (atom.resi < protein.helix[j][1])
                    continue;
                if (atom.resi > protein.helix[j][3])
                    continue;
                atom.ss = 'h';
                if (atom.resi == protein.helix[j][1])
                    atom.ssbegin = true;
                else if (atom.resi == protein.helix[j][3])
                    atom.ssend = true;
            }
        }

        return {
            'ok': true,
            'atoms': atoms,
            'bonds': bonds
        };
    };

    function parseAtom(line, atoms, atomSerialMap) {
        var x = parseFloat(line.substr(30, 7));
        var y = parseFloat(line.substr(38, 7));
        var z = parseFloat(line.substr(46, 7));
        var resi = parseInt(line.substr(22, 5));
        var chain = line.substr(21, 1);

        var e = trim(line.substr(76, 2)).toLowerCase(); //element ie 'C' for carbon

        if (e === '')
            e = trim(line.substr(12, 2)).toLowerCase();

        var atom = trim(line.substr(12, 2)).toLowerCase();

        var hflag;
        (line.substr(0, 6) === 'HETATM') ? hflag = true : hflag = false;

        var serial = line.substr(7, 4);
        atomSerialMap.put(serial, atoms.size);
        atoms.push({x: x, y: y, z: z, element: e, hflag: hflag, ss: 'c', atom: atom, resi: resi, chain: chain});
    }

    function parseConnect(line, bonds, bondHash, atomSerialMap) {
        var gap = 5; //the gap between the atom to be bonded
        var startIndex = 11; //the start of the first serial

        var startSerial = parseInt(line.substr(6, 5));
        var startAtom = atomSerialMap.get(startSerial);

        for (var i = 0; i < 4; ++i) {

            var lineIndex = startIndex + gap * i;
            var endSerial = parseInt(line.substr(lineIndex, gap));
            var endAtom = atomSerialMap.get(endSerial);

            if (endAtom) {

                var h = hash(startAtom, endAtom); // used to check for multiple bonds

                if (bondHash[ h ] === undefined) {

                    bondHash[ h ] = bonds.length;
                    bonds.push({from: startAtom - 1, to: endAtom - 1, count: 1});
                } else {

                    bonds[bondHash[h]].count += 1;
                }
            }
        }
    }

    function parseHelix(line, protein) {

        var startChain = line.substr(19, 1);
        var startResi = parseInt(line.substr(21, 4));
        var endChain = line.substr(31, 1);
        var endResi = parseInt(line.substr(33, 4));
        protein.helix.push([startChain, startResi, endChain, endResi]);
    }

    function trim(text) {

        return text.replace(/^\s\s*/, '').replace(/\s\s*$/, '');
    }

    function hash(s, e) {

        return 's' + Math.min(s, e) + 'e' + Math.max(s, e);
    }

    function defineCell(protein) {

        if (!protein.a) // check if undefined
            return;

        protein.ax = protein.a;
        protein.ay = 0;
        protein.az = 0;
        protein.bx = protein.b * Math.cos(Math.PI / 180.0 * protein.gamma);
        protein.by = protein.b * Math.sin(Math.PI / 180.0 * protein.gamma);
        protein.bz = 0;
        protein.cx = protein.c * Math.cos(Math.PI / 180.0 * protein.beta);
        protein.cy = protein.c * (Math.cos(Math.PI / 180.0 * protein.alpha) -
                Math.cos(Math.PI / 180.0 * protein.gamma)
                * Math.cos(Math.PI / 180.0 * protein.beta)
                / Math.sin(Math.PI / 180.0 * protein.gamma));
        protein.cz = Math.sqrt(protein.c * protein.c * Math.sin(Math.PI / 180.0 * protein.beta)
                * Math.sin(Math.PI / 180.0 * protein.beta) - protein.cy * protein.cy);
    }

    window.PDBParser = PDBParser;
})(window);
