var Erlang = require('plom-erlangify')
  , fs = require('fs');

//get some model components
var p = require('../noise/process.json')
  , l = require('../noise/link.json')
  , t = require('../noise/theta.json');

//Definition of the Erlang expansion
var def = [
  {from: 'I', to: 'I', rate: 'correct_rate(v)', shape: 2, rescale: 'v'}
];

//erlangify !    
var erlang = new Erlang(def);

fs.writeFileSync('process.json', JSON.stringify(erlang.ify(p)));
fs.writeFileSync('link.json', JSON.stringify(erlang.ify(l)));
fs.writeFileSync('theta.json', JSON.stringify(erlang.ify(t)));
