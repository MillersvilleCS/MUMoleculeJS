(function(window) {
    'use strict';

    var MoleculeGeometryBuilder = function( ) {

    };

    MoleculeGeometryBuilder.BONDS_LINES = 0;
    MoleculeGeometryBuilder.BONDS_CYLINDERS = 1;

    //MoleculeGeometryBuilder.ATOMS_CIRCLES = 0;
    MoleculeGeometryBuilder.ATOMS_SPHERES = 1;

    MoleculeGeometryBuilder.create = function(json, atomRenderType, bondRenderType, atomScale, bondThickness) {

        if (json.ok === false)
            throw new Exception("Cannon create model, invalid json");

        var model = new THREE.Object3D( );

        var atoms = json.atoms;
        var bonds = json.bonds;

        //createRibbons(model, atoms, 3, 5, 'p');
        //createRibbons(model, atoms, 3, 5, 'o');
        drawStrand(model, atoms, 6, 5, false, 3, 1.3, false, 3);

        switch (atomRenderType) {

            case MoleculeGeometryBuilder.ATOMS_SPHERES:
                createAtomsAsSpheres(atoms, atomScale, 16, model);
                break;
            default:
                createAtomsAsSpheres(atoms, atomScale, 16, model);
        }

        switch (bondRenderType) {

            case MoleculeGeometryBuilder.BONDS_LINES:
                createBondsAsLines(bonds, atoms, bondThickness, model);
                break;
            case MoleculeGeometryBuilder.BONDS_CYLINDERS:
                createBondsAsCylinders(bonds, atoms, bondThickness, model);
                break;
            default:
                createBondsAsLines(bonds, atoms, bondThickness, model);
        }

        return model;
    };

    function createRibbons(model, atoms, curveWidth, div, chainElement) {

        var points = [], colors = [];
        var currentChain, currentResi;

        for (var i = 0; i < atoms.length; i++) {

            var atom = atoms[i];

            if (atom.hflag)
                continue;

            if (atom.atom === chainElement) {

                if (currentChain !== atom.chain || currentResi + 1 !== atom.resi) {

                    createSmoothCurve(model, points, curveWidth, colors, div);
                    points = [];
                    colors = [];
                }

                points.push(new THREE.Vector3(atom.x, atom.y, atom.z));
                colors.push(getAtomColor(atom.element));
                currentChain = atom.chain;
                currentResi = atom.resi;
            }
        }

        createSmoothCurve(model, points, curveWidth, colors, div);
    }

    function createSmoothCurve(model, points, width, colors, div) {

        if (points.length === 0)
            return;

        var geo = new THREE.Geometry();
        var dividedPoints = subdivide(points, div);

        for (var i = 0; i < dividedPoints.length; i++) {

            var colorIndex = (i === 0) ? 0 : Math.round((i - 1) / div); //
            geo.vertices.push(dividedPoints[i]);
            geo.colors.push(colors[colorIndex]);
        }

        var lineMaterial = new THREE.LineBasicMaterial({linewidth: width});
        lineMaterial.vertexColors = true;
        var line = new THREE.Line(geo, lineMaterial);
        line.type = THREE.LineStrip;
        model.add(line);
    }

    function subdivide(_points, DIV) { // points as Vector3
        var ret = [];
        var points = _points;
        points = new Array(); // Smoothing test
        points.push(_points[0]);
        for (var i = 1, lim = _points.length - 1; i < lim; i++) {
            var p1 = _points[i], p2 = _points[i + 1];
            if (p1.smoothen)
                points.push(new TV3((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2));
            else
                points.push(p1);
        }
        points.push(_points[_points.length - 1]);

        for (var i = -1, size = points.length; i <= size - 3; i++) {
            var p0 = points[(i == -1) ? 0 : i];
            var p1 = points[i + 1], p2 = points[i + 2];
            var p3 = points[(i == size - 3) ? size - 1 : i + 3];
            var v0 = new TV3().sub(p2, p0).multiplyScalar(0.5);
            var v1 = new TV3().sub(p3, p1).multiplyScalar(0.5);
            for (var j = 0; j < DIV; j++) {
                var t = 1.0 / DIV * j;
                var x = p1.x + t * v0.x
                        + t * t * (-3 * p1.x + 3 * p2.x - 2 * v0.x - v1.x)
                        + t * t * t * (2 * p1.x - 2 * p2.x + v0.x + v1.x);
                var y = p1.y + t * v0.y
                        + t * t * (-3 * p1.y + 3 * p2.y - 2 * v0.y - v1.y)
                        + t * t * t * (2 * p1.y - 2 * p2.y + v0.y + v1.y);
                var z = p1.z + t * v0.z
                        + t * t * (-3 * p1.z + 3 * p2.z - 2 * v0.z - v1.z)
                        + t * t * t * (2 * p1.z - 2 * p2.z + v0.z + v1.z);
                ret.push(new TV3(x, y, z));
            }
        }
        ret.push(points[points.length - 1]);
        return ret;
    }

    function drawStrand(group, atoms, num, div, fill, coilWidth, helixSheetWidth, doNotSmoothen, thickness) {
        
        var currentChain, currentResi, currentCA;
        var prevCO = null, ss = null;
        
        var points = [];
        for (var k = 0; k < num; k++)
            points[k] = [];
        
        var colors = [];  

        for (var i = 0; i < atoms.length; i++) {

            var atom = atoms[i];

            if (atom.hflag)
                continue;

            if (atom.atom === 'ca') {

                if (currentChain !== atom.chain || currentResi + 1 !== atom.resi) {

                    for (var j = 0; !thickness && j < num; j++)
                        drawSmoothCurve(group, points[j], 1, colors, div);

                    if (fill)
                        drawStrip(group, points[0], points[num - 1], colors, div, thickness);

                    var points = [];

                    for (var k = 0; k < num; k++)
                        points[k] = [];

                    colors = [];
                    prevCO = null;
                    ss = null;
                }

                currentCA = new THREE.Vector3(atom.x, atom.y, atom.z);
                currentChain = atom.chain;
                currentResi = atom.resi;
                ss = atom.ss;
                colors.push(getAtomColor(atom.element));
                
            } else if(atom.atom === 'o') { 
                
                if (!currentCA) {
                    currentCA = new THREE.Vector3(atom.x, atom.y, atom.z);
                    continue;
                }
                
                var O = new THREE.Vector3(atom.x, atom.y, atom.z);
                O.sub(currentCA);
                O.normalize(); // can be omitted for performance
                O.multiplyScalar((ss == 'c') ? coilWidth : helixSheetWidth);

                if (prevCO != undefined && O.dot(prevCO) < 0)
                    O.negate();

                prevCO = O;

                for (var j = 0; j < num; j++) {
                    
                    var delta = -1 + 2 / (num - 1) * j;
                    var v = new THREE.Vector3(currentCA.x + prevCO.x * delta, currentCA.y + prevCO.y * delta, currentCA.z + prevCO.z * delta);

                    if (!doNotSmoothen && ss == 's')
                        v.smoothen = true;

                    points[j].push(v);
                }
            }
        }

        for (var j = 0; !thickness && j < num; j++)
            createSmoothCurve(group, points[j], 1, colors, div);

        if (fill)
            drawStrip(group, points[0], points[num - 1], colors, div, thickness);
    }

    function drawStrip(group, p1, p2, colors, div, thickness) {

        if ((p1.length) < 2)
            return;

        p1 = subdivide(p1, div);
        p2 = subdivide(p2, div);

        if (!thickness)
            return drawThinStrip(group, p1, p2, colors, div);

        var geo = new THREE.Geometry();
        var vs = geo.vertices, fs = geo.faces;
        var axis, p1v, p2v, a1v, a2v;

        for (var i = 0, lim = p1.length; i < lim; i++) {
            vs.push(p1v = p1[i]); // 0
            vs.push(p1v); // 1
            vs.push(p2v = p2[i]); // 2
            vs.push(p2v); // 3
            if (i < lim - 1) {
                var toNext = p1[i + 1].clone().sub(p1[i]);
                var toSide = p2[i].clone().sub(p1[i]);
                axis = toSide.crossSelf(toNext).normalize().multiplyScalar(thickness);
            }
            vs.push(a1v = p1[i].clone().add(axis)); // 4
            vs.push(a1v); // 5
            vs.push(a2v = p2[i].clone().add(axis)); // 6
            vs.push(a2v); // 7
        }
        var faces = [[0, 2, -6, -8], [-4, -2, 6, 4], [7, 3, -5, -1], [-3, -7, 1, 5]];
        for (var i = 1, lim = p1.length; i < lim; i++) {
            var offset = 8 * i, color = new THREE.Color(colors[Math.round((i - 1) / div)]);
            for (var j = 0; j < 4; j++) {
                var f = new THREE.Face4(offset + faces[j][0], offset + faces[j][1], offset + faces[j][2], offset + faces[j][3], undefined, color);
                fs.push(f);
            }
        }
        var vsize = vs.length - 8; // Cap
        for (var i = 0; i < 4; i++) {
            vs.push(vs[i * 2]);
            vs.push(vs[vsize + i * 2])
        }
        ;
        vsize += 8;
        fs.push(new THREE.Face4(vsize, vsize + 2, vsize + 6, vsize + 4, undefined, fs[0].color));
        fs.push(new THREE.Face4(vsize + 1, vsize + 5, vsize + 7, vsize + 3, undefined, fs[fs.length - 3].color));
        geo.computeFaceNormals();
        geo.computeVertexNormals(false);
        var material = new THREE.MeshLambertMaterial();
        material.vertexColors = THREE.FaceColors;
        var mesh = new THREE.Mesh(geo, material);
        mesh.doubleSided = true;
        group.add(mesh);
    }

    function drawThinStrip(group, p1, p2, colors, div) {
        var geo = new THREE.Geometry();
        for (var i = 0, lim = p1.length; i < lim; i++) {
            geo.vertices.push(p1[i]); // 2i
            geo.vertices.push(p2[i]); // 2i + 1
        }
        for (var i = 1, lim = p1.length; i < lim; i++) {
            var f = new THREE.Face4(2 * i, 2 * i + 1, 2 * i - 1, 2 * i - 2);
            f.color = new THREE.Color(colors[Math.round((i - 1) / div)]);
            geo.faces.push(f);
        }
        geo.computeFaceNormals();
        geo.computeVertexNormals(false);
        var material = new THREE.MeshLambertMaterial();
        material.vertexColors = THREE.FaceColors;
        var mesh = new THREE.Mesh(geo, material);
        mesh.doubleSided = true;
        group.add(mesh);
    }


    function createAtomsAsSpheres(atoms, scale, quality, model) {

        var geometry = new THREE.SphereGeometry(1, quality, quality);

        for (var i = 0; i < atoms.length; i++) {

            var atom = atoms[i];

            if (!atom.hflag) //hflag denotes a normal atom
                continue;

            var position = new THREE.Vector3(atom.x, atom.y, atom.z);

            var color = getAtomColor(atom.element);

            var radius = getAtomRadius(atom.element, scale);

            var material = new THREE.MeshLambertMaterial({color: color.getHex( )});

            var mesh = new THREE.Mesh(geometry.clone( ), material);

            mesh.scale.x *= radius;
            mesh.scale.y *= radius;
            mesh.scale.z *= radius;
            mesh.position = position;

            model.add(mesh);
        }
    }

    function createBondsAsLines(bonds, atoms, thickness, model) {

        for (var i = 0; i < bonds.length; i++) {

            var bond = bonds[i];

            var start = atoms[bond.from];
            var end = atoms[bond.to];
            var num = bond.count;

            var vertex1 = new THREE.Vector3(start.x, start.y, start.z);
            var vertex2 = new THREE.Vector3(end.x, end.y, end.z);

            var color1 = getAtomColor(start.element);
            var color2 = getAtomColor(end.element);

            var distVec = new THREE.Vector3(0, 0, 0);

            var vertex1Sub = vertex1.clone().sub(distVec);
            var vertex1Add = vertex1.clone().add(distVec);
            var vertex2Sub = vertex2.clone().sub(distVec);
            var vertex2Add = vertex2.clone().add(distVec);

            switch (num) {
                case 3:
                    model.add(GeometryUtil.createLine(vertex1Sub, vertex2Sub, color1, color2));
                    model.add(GeometryUtil.createLine(vertex1, vertex2, color1, color2));
                    model.add(GeometryUtil.createLine(vertex1Add, vertex2Add, color1, color2));
                    break;
                case 2:
                    model.add(GeometryUtil.createLine(vertex1Sub, vertex2Sub, color1, color2));
                    model.add(GeometryUtil.createLine(vertex1Add, vertex2Add, color1, color2));
                    break;
                default:
                    model.add(GeometryUtil.createLine(vertex1, vertex2, color1, color2));
            }
        }
    }

    function createBondsAsCylinders(bonds, atoms, distanceApart, model) {

        for (var i = 0; i < bonds.length; i++) {

            var bond = bonds[i];

            var start = atoms[bond.from];
            var end = atoms[bond.to];
            var num = bond.count;

            var vertex1 = new THREE.Vector3(start.x, start.y, start.z);
            var vertex2 = new THREE.Vector3(end.x, end.y, end.z);

            var color1 = getAtomColor(start.element);
            var color2 = getAtomColor(end.element);

            var distVec = new THREE.Vector3(distanceApart, distanceApart, distanceApart);

            var vertex1Sub = vertex1.clone().sub(distVec);
            var vertex1Add = vertex1.clone().add(distVec);
            var vertex2Sub = vertex2.clone().sub(distVec);
            var vertex2Add = vertex2.clone().add(distVec);

            switch (num) {
                case 3:
                    model.add(GeometryUtil.createCylinder(vertex1Sub, vertex2Sub, color1, color2));
                    model.add(GeometryUtil.createCylinder(vertex1, vertex2, color1, color2));
                    model.add(GeometryUtil.createCylinder(vertex1Add, vertex2Add, color1, color2));
                    break;
                case 2:
                    model.add(GeometryUtil.createCylinder(vertex1Sub, vertex2Sub, color1, color2));
                    model.add(GeometryUtil.createCylinder(vertex1Add, vertex2Add, color1, color2));
                    break;
                default:
                    model.add(GeometryUtil.createCylinder(vertex1, vertex2, color1, color2));
            }
        }
    }

    function getAtomColor(element) {

        var atomColor = atomColors[element.toString()];

        var r = atomColor[0] / 255;
        var g = atomColor[1] / 255;
        var b = atomColor[2] / 255;

        var color = new THREE.Color();

        return color.setRGB(r, g, b);
    }

    function getAtomRadius(element, scale) {

        var radius = atomRadii[element.toString()];

        if (radius === undefined)
            radius = atomRadii["default"];

        radius *= scale;

        return radius;
    }

    var atomColors = {
        'h': [255, 255, 255],
        'he': [217, 255, 255],
        'li': [204, 128, 255],
        'be': [194, 255, 0],
        'b': [255, 181, 181],
        'c': [144, 144, 144],
        'n': [48, 80, 248],
        'o': [255, 13, 13],
        'f': [144, 224, 80],
        'ne': [179, 227, 245],
        'na': [171, 92, 242],
        'mg': [138, 255, 0],
        'al': [191, 166, 166],
        'si': [240, 200, 160],
        'p': [255, 128, 0],
        's': [255, 255, 48],
        'cl': [31, 240, 31],
        'ar': [128, 209, 227],
        'k': [143, 64, 212],
        'ca': [61, 255, 0],
        'sc': [230, 230, 230],
        'ti': [191, 194, 199],
        'v': [166, 166, 171],
        'cr': [138, 153, 199],
        'mn': [156, 122, 199],
        'fe': [224, 102, 51],
        'co': [240, 144, 160],
        'ni': [80, 208, 80],
        'cu': [200, 128, 51],
        'zn': [125, 128, 176],
        'ga': [194, 143, 143],
        'ge': [102, 143, 143],
        'as': [189, 128, 227],
        'se': [255, 161, 0],
        'br': [166, 41, 41],
        'kr': [92, 184, 209],
        'rb': [112, 46, 176],
        'sr': [0, 255, 0],
        'y': [148, 255, 255],
        'zr': [148, 224, 224],
        'nb': [115, 194, 201],
        'mo': [84, 181, 181],
        'tc': [59, 158, 158],
        'ru': [36, 143, 143],
        'rh': [10, 125, 140],
        'pd': [0, 105, 133],
        'ag': [192, 192, 192],
        'cd': [255, 217, 143],
        'in': [166, 117, 115],
        'sn': [102, 128, 128],
        'sb': [158, 99, 181],
        'te': [212, 122, 0],
        'i': [148, 0, 148],
        'xe': [66, 158, 176],
        'cs': [87, 23, 143],
        'ba': [0, 201, 0],
        'la': [112, 212, 255],
        'ce': [255, 255, 199],
        'pr': [217, 255, 199],
        'nd': [199, 255, 199],
        'pm': [163, 255, 199],
        'sm': [143, 255, 199],
        'eu': [97, 255, 199],
        'gd': [69, 255, 199],
        'tb': [48, 255, 199],
        'dy': [31, 255, 199],
        'ho': [0, 255, 156],
        'er': [0, 230, 117],
        'tm': [0, 212, 82],
        'yb': [0, 191, 56],
        'lu': [0, 171, 36],
        'hf': [77, 194, 255],
        'ta': [77, 166, 255],
        'w': [33, 148, 214],
        're': [38, 125, 171],
        'os': [38, 102, 150],
        'ir': [23, 84, 135],
        'pt': [208, 208, 224],
        'au': [255, 209, 35],
        'hg': [184, 184, 208],
        'tl': [166, 84, 77],
        'pb': [87, 89, 97],
        'bi': [158, 79, 181],
        'po': [171, 92, 0],
        'at': [117, 79, 69],
        'rn': [66, 130, 150],
        'fr': [66, 0, 102],
        'ra': [0, 125, 0],
        'ac': [112, 171, 250],
        'th': [0, 186, 255],
        'pa': [0, 161, 255],
        'u': [0, 143, 255],
        'np': [0, 128, 255],
        'pu': [0, 107, 255],
        'am': [84, 92, 242],
        'cm': [120, 92, 227],
        'bk': [138, 79, 227],
        'cf': [161, 54, 212],
        'es': [179, 31, 212],
        'fm': [179, 31, 186],
        'md': [179, 13, 166],
        'no': [189, 13, 135],
        'lr': [199, 0, 102],
        'rf': [204, 0, 89],
        'db': [209, 0, 79],
        'sg': [217, 0, 69],
        'bh': [224, 0, 56],
        'hs': [230, 0, 46],
        'mt': [235, 0, 38],
        'ds': [235, 0, 38],
        'rg': [235, 0, 38],
        'cn': [235, 0, 38],
        'uut': [235, 0, 38],
        'uuq': [235, 0, 38],
        'uup': [235, 0, 38],
        'uuh': [235, 0, 38],
        'uus': [235, 0, 38],
        'uuo': [235, 0, 38]
    };

    var atomRadii = {
        'h': 1.2,
        'li': 1.82,
        'na': 2.27,
        'k': 2.75,
        'c': 1.7,
        'n': 1.55,
        'o': 1.52,
        'f': 1.47,
        'p': 1.80,
        's': 1.80,
        'cl': 1.75,
        'br': 1.85,
        'se': 1.90,
        'zn': 1.39,
        'cu': 1.4,
        'ni': 1.63,
        'default': 1.5
    };

    window.MoleculeGeometryBuilder = MoleculeGeometryBuilder;
})(window);
