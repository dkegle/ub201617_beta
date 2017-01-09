# Evolucija virusov
Skupina ub201617_beta

# Uvod

Cilj naloge je bil identificirati gene, ki vsebujejo informacije o tem kako se virus razvija med izbruhom,
opraviti globalno poravnano (global alignment) na izbranih genih in izračunati matriko razdalj med vzorci (edit distance).
S pomočjo znanih algoritmov gručenja (Hierarchical clustering, Neighbour joining, Multidimensional scaling ...) 
je bilo potrebno poročati o zanimivih ugotovitvah. Enega od teh algoritmov je bilo potrebno tudi implementirati.

# Metode in Rezultati

## Izračun matrike razdalje urejanja (edit distance)

Najprej smo identificirali gene in njihove velikosti. Ugotovili smo, da virus Ebole vsebuje 7 genov. Od tega so štirje geni krajši (VP24, VP30, VP35 in VP40), dva srednje velika (GP in NP) in eden daljši (L). Krajši štirje so dolgi približno 1500 bp, srednji približno 2500 do 3000 bp in najdaljši približno 7000 bp. 

Nato smo izbrali dva naključna vzorca in na vsakem genu z uporabo dinamičnega programiranja izvedli globalno poravnano (global alignment). Ker nas je zanimala samo razdalja urejanja (edit distance), traceback-a nismo potrebovali. Ugotovili smo, da za izračun poravnave enega para genov potrebujemo nekaj minut. Izračuna za najdaljši gen L nismo mogli izvesti na računalniku z 8GB spomina, saj je tabela dinamičnega programiranja postala prevelika.

Celotna matrika razdalj, ki se je uporabljala za gručenje, je vsota matrik razdalj za posamezne gene. Zaradi omejitve spomina smo gen L izpustili. Za pohitritev izračunov smo uporabili vzporedno računanje na več procesorskih jedrih. Za izračun celotne matrike razdalj je bilo potrebno približno 20 ur.

## Hierarhično gručenje

Za ocenjevanje podobnosti med različnimi genomi virusa ebole smo uporabili hierarhično gručenje.
Pri tem smo razdaljo med dvema genomoma definirali kot vsoto ocen poravnav (alignment score) med
5 geni iz vsakega genoma, za katere smo izračunali oceno poravnave.
Pri računanju razdalje med gručami pa smo uporabili povprečno razdaljo, tj. razdalja med dvema gručama je
povprečje razdalj med pari primerkov iz obeh gruč.

Postopek hierarhičnega gručenja smo sprogramirali sami, za izris dendrograma pa uporabili funkcijo iz
paketa scipy.

<!--![](hierarchial_dendrogram.png)-->
<a href="https://raw.githubusercontent.com/dkegle/ub201617_beta/master/hierarchial_dendrogram.png">
<img border="0" src="hierarchial_dendrogram.png" width="1024" height="768">
</a>

Na dendrogramu za vsak genom označimo državo, v kateri se je pojavil. Tako vidimo, da so v zgornji (redeči) skupini virusi opazno različni od
ostalih, pojavljajo pa se predvsem v Gvineji.
Opazimo še skupino zelo podobnih virusov, ki se pojavljajo predvsem v Sierri Leone, pa tudi v Gvineji in Liberiji. Sklepamo,
da so ti virusi v Gvinejo in Liberijo prišli iz Sierre Leone.
Zanimivo je tudi, da virusi iz Demokratične republike Kongo niso podobni nobenemu virusu iz sosednjih držav.

## MDS

Poleg dendrograma, smo za vizualizacijo izračunanih razdalj uporabili tudi MDS (Večrazsežnostno lestvičenje).
MDS podatke o podobnosti/različnosti preslika v razdalje med točkami v geometričnem prostoru. Da je vizualizacija kar se
da jasna, smo razdalje med genomi preslikali v 2-dimenzionalni prostor. Pozicije (koordinate) genomov smo izračunali z
uporabo funkcije mainfold.MDS iz knjižnjice scikit-learn. Začetne pozicije točk smo izbrali naključno, nato pa razdalje
med njimi optimizirali z metričnim MDS.

Dobljene točke smo prikazali na scatterplotu (razsevni diagram) in jih obarvali glede na državo pojavitve.
Barve smo izbrali tako, da so države, ki so si geografsko blizu, podobno obarvane.

![](mds_full.png)

Opazimo tesno gručo genomov iz Gabona in Demokratične republike Kongo. 

![](mds_zoom2.png)

V bližini tudi skupino zelo podobnih genomov
iz Sierre Leone, v bližini katere se nahajajo tudi vzorci iz bližnjih Liberie in Gvineje.
Opazno različna od vseh prej omenjenih skupin pa je gruča genomov, pretežno iz Gvineje, ki je obarvana rdeče in izrisana
na povsem drugi strani diagrama.
Med obravnavanimi genomi je tudi nekaj vzorcev iz Nigerije, ki pa glede na MDS ne spadajo v nobeno od identificiranih
gruč - verjetno tudi zaradi geografske oddaljenosti.

## Neighbour joining

Neighbour joining je metoda gručenja za izgradnjo genetskih dreves. Metodo sta ustvarila Naruya Saitou in Masatoshi Nei leta 1987.
Običajno se metoda uporablja na drevesih, ki temeljijo na DNA ali proteinskih podatkih. Metoda zahteva razdalje med pari s katerimi lahko zgradimo drevo.
Končno drevo nam prikaže kateri pari so si med seboj najbolj dobobni. 

![](neighbour_joining.png)

Na zgornjem drevesu vidimo, da so primerki iz enakih držav praviloma blizu med seboj. Naprimer vsi trije primerki iz Gabona so si med seboj zelo blizu.
Enako velja tudi za oba iz Nigerije. Iz Demokratične republike Congo so 4je primerki zelo skupaj in so si med seboj podobni s primerki iz Gabona, eden pa malo izstopa
in je na grafu daleč stran. Najbolj razpršeni so primerki iz Liberije in Sierre Leone, kar je po svoje razumljivo, saj jih je največ.


