; Last Updated: 6-11-20
; Creator: David Elzinga

; We need the R extension in order to use beta distributions,
; so enable the R extension... we also need the csv extension to
; read in our csv
extensions[ csv ]

; We have two types of pollen, Columbia (blue) and Landsberg (red)
; but they are grouped under one breed, pollen, with color properties.
breed[pollen pollens]

; Globals.
globals[stigma-length num-pollen new-prop-col new-prop-ler first-germ-flag attract-mult LHS-params-table Betas-table params r3C r6C r9C r24C r3L r6L r9L r24L]

; Agent variables.
pollen-own[ germ-check energy effective-length parent ]

; Patch variables.
patches-own[ starch attraction repelling ]

;-------- Set up --------;
to setup

  ; Basic clean up
  clear-all
  reset-ticks

  if run-type = "LHS" or run-type = "best_mean_params" or run-type = "best_distr_params"
  [
  ; Input LHS parameters.
    if run-type = "LHS" [
      set LHS-params-table csv:from-file "LHS.csv"]
    if run-type = "best_mean_params" [
      set LHS-params-table csv:from-file "best_mean_params.csv"]
    if run-type = "best_distr_params" [
      set LHS-params-table csv:from-file "best_distr_params.csv"]

  ; Assign our parameters according to the LHS. Remember, NetLogo indexing starts at zero! STARTS AT ZERO
  set params item LHS-paramcombo-counter LHS-params-table
  set chem-atr-rad-ler round(item 0 params)
  set chem-atr-rad-col round(item 1 params)
  set initial-starch item 2 params
  set init-energy-ler item 3 params
  set init-energy-col item 4 params
  set max-starch-harvest-col item 5 params
  set max-starch-harvest-ler item 6 params

  (ifelse item 7 params <= 0.3333
    [set replenish-style "no-replenish"]
    item 7 params <= 0.6666
    [set replenish-style "replenish-all"]
    item 7 params < 1
    [set replenish-style "replenish-empty"])

  set fert-thresh-col item 8 params
  set fert-thresh-ler item 9 params
  set prop-fert-cutoff item 10 params
  set freq-starch-replenish round(item 11 params)
  set prop-starch-replenish item 12 params
  set step-size item 13 params
  set movement-cost-ler item 14 params
  set movement-cost-col item 15 params
  set H1P1 item 16 params
  ]

  ; Define the length of the style and stigma. Entire pistil is about 3.5mm long, stigma is about 0.37mm long, the style is about 0.14mm long,
  ; and the transmitting trac is about 0.17mm wide. The world is 146 patches tall (146 at the top of the world). To convert from mm to patches, then, 146/3.5 * mm
  let style-length round(max-pycor / 3.5 * 0.14)
  set stigma-length round(max-pycor / 3.5 * 0.37)
  let trans-trac-width round(max-pycor / 3.5 * 0.17)
  let ovary-length max-pycor - style-length - stigma-length

  ; Create the ovules. The user specifies the total number of ovules (an even number)
  ; that are distributed equally on the the left and right side of the
  ; pistil. The number of ovules comes from Swanson's study. The spacing
  ; of the ovules is such that 60 can fit in the ovary (30 on each side).
  let num-ovules 60
  let ovule-spacing floor (ovary-length / (num-ovules / 2))
  let max-chem-rad max list chem-atr-rad-ler chem-atr-rad-col
  let counter 0
  while [counter < num-ovules / 2]
  [
    ; Ovules begin immedately following the stigma (at the top of the style). (separated so that 60 ovules can fit)
    ask patches with [pycor = max-pycor - stigma-length - style-length - (ovule-spacing * counter) and ((pxcor > max-pxcor - 1) or (pxcor < min-pxcor + 1))]
          [set pcolor yellow

           ; Patches within the chemoattractant range of at least one of the colors
           ; should be marked green. These patches are not the yellow patch.
           ask patches in-radius max-chem-rad with [pcolor != yellow]
                  [set pcolor green]
           ]
    set counter counter + 1
  ]

  ; Define the total amount of pollen grains to be applied to the stigma
  ; which comes from another Swanson proposal
  set num-pollen random-normal 1040 236

  ; Create the pollen! First make the red (landsberg) pollen, then blue (col).
  ; If we are running just a single accession, adjust the proportions.

  if accession = "ler-only"
  [set prop-ler 1]
  if accession = "col-only"
  [set prop-ler 0]

  create-pollen num-pollen * prop-ler
  [
    set color red
    set parent 0
    set germ-check 0
    ifelse random-float 1 < 0.5
    [setxy ((random ((trans-trac-width + 1) / 2)) + max-pxcor / 2) ((random stigma-length + 1) + ovary-length + style-length)]
    [setxy ((random (- (trans-trac-width + 1 )/ 2)) + max-pxcor / 2) ((random stigma-length + 1) + ovary-length + style-length)]
  ]

  create-pollen num-pollen - count pollen with [color = red]
  [
    set color blue
    set parent 0
    set germ-check 0
    ifelse random-float 1 < 0.5
    [setxy ((random ((trans-trac-width + 1) / 2)) + max-pxcor / 2) ((random stigma-length + 1) + ovary-length + style-length)]
    [setxy ((random (- (trans-trac-width + 1 )/ 2)) + max-pxcor / 2) ((random stigma-length + 1) + ovary-length + style-length)]
  ]

  set-default-shape pollen "circle"
  ; *** = where do these dimensions come from?
  ; why are the ovulues spaced four apart?, why is there a 0.5 shift in pollen placement?


  ; Now that our pollen is created, we kill off the pollen that won't germinate (remember,
  ; this encompases both pollen that aren't viable, and won't germinate).
  ;set Betas-table csv:from-file "Betas.csv"
  ;let beta-choice random 1000000
  ;let beta-row item beta-choice Betas-table

  ;let col-germ item 0 beta-row
  ;let ler-germ item 1 beta-row

  ;ask pollen
  ;[
  ;  if color = red and random-float 1 > ler-germ
  ;  [ die ]

  ;  if color = blue and random-float 1 > col-germ
  ;  [ die ]
  ;]

  ; All patches that are not yellow have some initial value of starch,
  ; assign that value now.
  ask patches with [pcolor != yellow]
  [set starch initial-starch]

end


;-------- go --------;

to go

  ; To accurately measure time in the experiment, once the first pollen tube
  ; germinates the ticks are reset.
  if (count pollen with [germ-check = 1]) > 0 and first-germ-flag = 0
  [
    set first-germ-flag 1
  ]

  ; Ask patches to replenish starch (if necessary)
  if replenish-style = "replenish-all" or replenish-style = "replenish-empty"
  [ ask patches
    [ replenish-starch ]
  ]

  ; Germinate some pollen.
  germ

  ; Update the attractant multiplier
  update-attract-mult

  ; Ask the pollen that have germinated to harvet, move, and fertilize.
  ask pollen with [germ-check = 1 and parent = 0 and ycor > 0]
  [
    harvest
    move
    fertilize
  ]

  ; Record that data!
  if first-germ-flag = 1 and (ticks = 180 or ticks = 360 or ticks = 540 or ticks = 1440)
  [record-data]

  ; Advance time (if we haven't reached 1440mins = 24hrs)!
  ifelse ticks < 1440
  [tick]
  [stop]
end

;-------- subroutines --------;
to germ

  ; First, germinate some columbia pollen. This quantity is based off the derivative
  ; of the graph that describes the change in proportion of pollen that have germinated
  ; (which we dictate must be at least zero), that is,

  let col-germ-table csv:from-file "col_germ.csv"
  let ler-germ-table csv:from-file "ler_germ.csv"


  ifelse ticks <= 40
  [
     set new-prop-col item 0 (item (ticks) col-germ-table)
     set new-prop-ler item 0 (item (ticks) ler-germ-table)
  ]
  [
    set new-prop-col 0
    set new-prop-ler 0
  ]

  let dif-col max list ((new-prop-col * (count pollen with [color = blue])) - (count pollen with [color = blue and germ-check = 1])) 0
  let dif-ler max list ((new-prop-ler * (count pollen with [color = red])) - (count pollen with [color = red and germ-check = 1])) 0

  show dif-col
  show dif-ler

  ask up-to-n-of dif-col (pollen with [color = blue and germ-check = 0])
  [
    set germ-check 1
    set energy init-energy-col
  ]

  ask up-to-n-of dif-ler (pollen with [color = red and germ-check = 0])
  [
    set germ-check 1
    set energy init-energy-col
  ]
end


to harvest

  ; Based on their color, the pollen should harvest up to the max
  ; amount of starch they can, and the patch should have its energy reduced
  ; in conjunction.
  if starch > 0
  [
    ifelse color = blue ; it's col!
    [
      let harvest-amount min list starch max-starch-harvest-col
      set energy energy + harvest-amount
      set starch starch - harvest-amount
    ]

    ; it's ler!
    [
      let harvest-amount min list starch max-starch-harvest-ler
      set energy energy + harvest-amount
      set starch starch - harvest-amount
    ]

  ]
end


to move

  let chem-atr-rad 0
  let movement-cost 0
  let move-flag 0
  let pollen-xcor xcor
  let pollen-ycor ycor
  ifelse color = blue
  [
    set chem-atr-rad chem-atr-rad-col
    set movement-cost movement-cost-col
  ]
  [
    set chem-atr-rad chem-atr-rad-ler
    set movement-cost movement-cost-ler
  ]

  ; Begin by facing a random patch looking down
  facexy  random-xcor 0

  ; Suppose that the pollen tube is not on an unfertilized ovule (otherwise it
  ; shouldn't move), and  there is an unfertilized ovule within the
  ; chemoattractant range (which depends upon the species of pollen).
  (ifelse count patches with [pcolor = yellow] in-radius chem-atr-rad >= 1 and count patches with [pcolor = yellow] in-radius 0 = 0
  [

   ; First, determine the potential patches that are attracting this pollen.
   let attracting-patches patch-set patches with [pcolor = yellow] in-radius chem-atr-rad

   ; Now we assign a certain attraction to each of the attracting-patches.
    ask attracting-patches
    [
     ; Define the attraction to each attracting patch, it should be in [0,1]
      set attraction attract-mult * (chem-atr-rad - (distancexy pollen-xcor pollen-ycor)) / chem-atr-rad
    ]

    ; Face one of the attracting-patches with the maximum attraction value.
    let most-attract-patch max-one-of attracting-patches [attraction]

    ; With probability of attraction, we face the most attractive patch. Otherwise,
    ; we'll just face a random patch downward.
    if [attraction] of most-attract-patch >= random-float 1
    [ face most-attract-patch ]
  ]

  ; Suppose there aren't any unfertilized ovules. If there are fertilized ovules,
  ; we need to repel them. We assume the chemorepellent range to be similar
  ; to the chemoattractant range.
  count patches with [pcolor = blue or pcolor = red] in-radius chem-atr-rad >= 1
  [

    ; First, determine the potential patches that are repelling this pollen.
    let repelling-patches patch-set patches with [pcolor = blue or pcolor = red] in-radius chem-atr-rad

    ; Now we assign a certain repllent to each of the attracting-patches.
    ask repelling-patches
    [
      ; Define the repellent to each repelling patch, it should be in [0,1]
       set repelling attract-mult * (chem-atr-rad - (distancexy pollen-xcor pollen-ycor)) / chem-atr-rad
    ]

    ; Face one of the repelling-patches with the maximum repelling value.
    let most-repelling-patch max-one-of repelling-patches [repelling]

    ; With probability of repelling, we face away (in x) from the most repelling patch. Otherwise,
    ; we'll just face a random patch downward.
    if [repelling] of most-repelling-patch >= random-float 1
     [ facexy  (max-pxcor - [pxcor] of most-repelling-patch) 0 ]
  ]

  ; Suppose there are not any attracting patches (besides ovules), and no repelling
  ; forces. Then, if they are on an ovlue, they shall not move. Put a movement flag of
  ; 1 on.
  count patches with [pcolor = yellow] in-radius 0 > 0
  [
      set move-flag 1
  ]
  ; If anything else happens, the pollen should just face a random direction and go there.
  [
      facexy  random-xcor 0
  ])

  ; Now let's move all pollen that don't have a movement flag of 1 and have sufficient energy
  if move-flag = 0 and energy > movement-cost
  [
   fd step-size
   set energy energy - movement-cost
  ]
end

to fertilize

  ; Check to see if we have a col with enough energy to fertilize at the right spot.
  if pcolor = yellow and (color = blue and energy >= fert-thresh-col)
  [
    ; If there is fertilization, we need to remove the ovule.
    set pcolor blue
    ; This is just for us to easily see that we have a blue ovule.
    ask neighbors
      [set pcolor blue]


    ; Now we have to update the chemoattractants to become chemorepellents.
    ; Since we are assuming that chemoattractants and chemorepellants don't overlap,
    ; then we just have to change the color of the chemoattracants from green.
    ; The actual implementation is a little more foregiving than that here, we
    ; ask all the attractant regions that aren't attractant regions for other
    ; ovules to turn to chemorepllent regions.
    let max-chem-rad max list chem-atr-rad-ler chem-atr-rad-col
    ask patches with [pcolor = green] in-radius (max-chem-rad + 1)
    [
      if count patches with [pcolor = yellow] in-radius max-chem-rad = 0
      [set pcolor gray]
    ]

    set parent 1
  ]

  ; Now we repeat the process for ler.
  ; Check to see if we have a col with enough energy to fertilize at the right spot.
  if pcolor = yellow and (color = red and energy >= fert-thresh-ler)
  [
    ; If there is fertilization, we need to remove the ovule.
    set pcolor red
    ; This is just for us to easily see that we have a red ovule.
    ask neighbors
      [set pcolor red]


    ; Now we have to update the chemoattractants to become chemorepellents.
    ; Since we are assuming that chemoattractants and chemorepellants don't overlap,
    ; then we just have to change the color of the chemoattracants from green.
    ; The actual implementation is a little more foregiving than that here, we
    ; ask all the attractant regions that aren't attractant regions for other
    ; ovules to turn to chemorepllent regions.
    let max-chem-rad max list chem-atr-rad-ler chem-atr-rad-col
    ask patches with [pcolor = green] in-radius (max-chem-rad + 1)
    [
      if count patches with [pcolor = yellow] in-radius max-chem-rad = 0
      [set pcolor gray]
    ]

    set parent 1
  ]
end

to replenish-starch

  ; Due to diminishing returns, there is likely some threshold at which
  ; producing starch is not worth it for the female. Only perform the following
  ; actions if we have not met that threshold yet. Also, if we are going to replenish
  ; starch, it can only happen at a certain frequency in ticks.
  let num-ovules 60
  if ticks mod freq-starch-replenish = 0 and (num-ovules - count patches with [pcolor = yellow]) / num-ovules <= prop-fert-cutoff
  [
    ; There are two circumstances underwhich we could have reached this point.
    ; Either, (1) we desire to replenish only the empty patches, or (2) we
    ; desire to replenish all patches.

    ifelse replenish-style = "replenish-all"
    [
      set starch min list initial-starch (starch + prop-starch-replenish * initial-starch)
    ]

    ; if we end up in this block, we need to only refill patches that are empty
    [
      if starch = 0
      [
        set starch min list initial-starch (starch + prop-starch-replenish * initial-starch)
      ]
    ]
  ]

end

to update-attract-mult

  ifelse chem-attrac-sat-fcn = "Holling-Type-1" and first-germ-flag = 1
  [
    set attract-mult min list (ticks * H1P1) 1
  ]
  ; If we haven't reached the initial germination, there's no attractants.
  [
    set attract-mult 0
  ]

end

to record-data

  ; Update the effective length of the pollen tubes.
  ask pollen
  [
    ifelse ycor > (max-pycor - stigma-length)
    [set effective-length 0]
    [set effective-length ((max-pycor - stigma-length) - ycor)]
  ]

  (ifelse ticks = 180
  [
      set r3C [(effective-length) *  3.5 / max-pycor] of pollen with [color = blue]
      set r3L [(effective-length) *  3.5 / max-pycor] of pollen with [color = red]
  ]
  ticks = 360
  [
      set r6C [(effective-length) *  3.5 / max-pycor] of pollen with [color = blue]
      set r6L [(effective-length) *  3.5 / max-pycor] of pollen with [color = red]
  ]
  ticks = 540
  [
      set r9C [(effective-length) *  3.5 / max-pycor] of pollen with [color = blue]
      set r9L [(effective-length) *  3.5 / max-pycor] of pollen with [color = red]
  ]
  ticks = 1440
  [
      set r24C [(effective-length) *  3.5 / max-pycor] of pollen with [color = blue]
      set r24L [(effective-length) *  3.5 / max-pycor] of pollen with [color = red]
  ])

end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
278
607
-1
-1
4.0
1
10
1
1
1
0
0
0
1
0
14
0
146
0
0
1
ticks
30.0

SLIDER
4
81
176
114
chem-atr-rad-ler
chem-atr-rad-ler
0
6
2.0
1
1
NIL
HORIZONTAL

SLIDER
4
116
176
149
chem-atr-rad-col
chem-atr-rad-col
0
6
6.0
1
1
NIL
HORIZONTAL

SLIDER
4
46
176
79
prop-ler
prop-ler
0
1
0.0
0.01
1
NIL
HORIZONTAL

BUTTON
5
10
68
43
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
4
151
177
184
initial-starch
initial-starch
0
100
31.6899522112288
1
1
NIL
HORIZONTAL

BUTTON
82
10
145
43
go 
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
4
186
177
219
init-energy-ler
init-energy-ler
0
100
67.4958183313315
1
1
NIL
HORIZONTAL

SLIDER
4
221
177
254
init-energy-col
init-energy-col
0
100
88.8139102187972
1
1
NIL
HORIZONTAL

SLIDER
4
256
177
289
max-starch-harvest-col
max-starch-harvest-col
0
20
7.69548030113114
0.5
1
NIL
HORIZONTAL

SLIDER
4
291
177
324
max-starch-harvest-ler
max-starch-harvest-ler
0
20
1.71128080777732
0.5
1
NIL
HORIZONTAL

CHOOSER
5
328
179
373
replenish-style
replenish-style
"no-replenish" "replenish-all" "replenish-empty"
2

SLIDER
5
375
180
408
fert-thresh-col
fert-thresh-col
0
20
18.795037799338104
1
1
NIL
HORIZONTAL

SLIDER
5
409
181
442
fert-thresh-ler
fert-thresh-ler
0
20
7.227915370651631
1
1
NIL
HORIZONTAL

SLIDER
5
444
181
477
prop-fert-cutoff
prop-fert-cutoff
0
1
0.48494982383551705
0.01
1
NIL
HORIZONTAL

SLIDER
6
479
181
512
freq-starch-replenish
freq-starch-replenish
0
720
450.0
1
1
NIL
HORIZONTAL

SLIDER
6
513
182
546
prop-starch-replenish
prop-starch-replenish
0
1
0.5420918952587771
0.01
1
NIL
HORIZONTAL

CHOOSER
293
14
431
59
chem-attrac-sat-fcn
chem-attrac-sat-fcn
"Holling-Type-1"
0

SLIDER
443
13
581
46
H1P1
H1P1
0
10
5.67960851642119
0.01
1
NIL
HORIZONTAL

SLIDER
293
62
465
95
step-size
step-size
0
5
0.27545987671797806
0.01
1
NIL
HORIZONTAL

SLIDER
293
100
465
133
movement-cost-ler
movement-cost-ler
0
10
7.33829335236743
1
1
NIL
HORIZONTAL

SLIDER
294
139
466
172
movement-cost-col
movement-cost-col
0
10
1.6251412069935602
1
1
NIL
HORIZONTAL

SLIDER
297
261
503
294
LHS-paramcombo-counter
LHS-paramcombo-counter
0
10000
91.0
1
1
NIL
HORIZONTAL

CHOOSER
297
216
454
261
run-type
run-type
"LHS" "manual" "best_distr_params" "best_mean_params"
3

CHOOSER
299
303
437
348
accession
accession
"both" "ler-only" "col-only"
0

BUTTON
707
147
770
180
go
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
311
379
526
424
# of Unfert. Ovules
count patches with [pcolor = yellow]
17
1
11

MONITOR
316
439
416
484
# of Ler Ovules
count patches with [pcolor = red] / 6
17
1
11

MONITOR
316
493
416
538
# of Col Ovules
count patches with [pcolor = blue] / 6
17
1
11

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Calibrate_LHS_Single_Accession" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1441"/>
    <metric>r3C</metric>
    <metric>r6C</metric>
    <metric>r9C</metric>
    <metric>r24C</metric>
    <metric>r3L</metric>
    <metric>r6L</metric>
    <metric>r9L</metric>
    <metric>r24L</metric>
    <metric>count patches with [pcolor = yellow]</metric>
    <metric>count patches with [pcolor = blue] / 6</metric>
    <metric>count patches with [pcolor = red] / 6</metric>
    <enumeratedValueSet variable="step-size">
      <value value="0.66"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-ler">
      <value value="18.975313187193"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-col">
      <value value="3.46683433710825"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accession">
      <value value="&quot;ler-only&quot;"/>
      <value value="&quot;col-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-type">
      <value value="&quot;LHS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-attrac-sat-fcn">
      <value value="&quot;Holling-Type-1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-starch">
      <value value="92.6892106482177"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-col">
      <value value="2.82183789445632"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-col">
      <value value="0.757218941401608"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-col">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="H1P1">
      <value value="5.34596861598446"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-ler">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-ler">
      <value value="5.64709511672291"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="LHS-paramcombo-counter" first="0" step="1" last="9999"/>
    <enumeratedValueSet variable="prop-ler">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replenish-style">
      <value value="&quot;replenish-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-ler">
      <value value="89.1021788475929"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-ler">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-fert-cutoff">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-col">
      <value value="9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Calibrate_LHS_Double_Accession" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1441"/>
    <metric>r3C</metric>
    <metric>r6C</metric>
    <metric>r9C</metric>
    <metric>r24C</metric>
    <metric>r3L</metric>
    <metric>r6L</metric>
    <metric>r9L</metric>
    <metric>r24L</metric>
    <metric>count patches with [pcolor = yellow]</metric>
    <metric>count patches with [pcolor = blue] / 6</metric>
    <metric>count patches with [pcolor = red] / 6</metric>
    <enumeratedValueSet variable="step-size">
      <value value="0.66"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-ler">
      <value value="18.975313187193"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-col">
      <value value="3.46683433710825"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accession">
      <value value="&quot;both&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-type">
      <value value="&quot;LHS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-attrac-sat-fcn">
      <value value="&quot;Holling-Type-1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-starch">
      <value value="92.6892106482177"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-col">
      <value value="2.82183789445632"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-col">
      <value value="0.757218941401608"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-col">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="H1P1">
      <value value="5.34596861598446"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-ler">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-ler">
      <value value="5.64709511672291"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="LHS-paramcombo-counter" first="0" step="1" last="9999"/>
    <enumeratedValueSet variable="prop-ler">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replenish-style">
      <value value="&quot;replenish-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-ler">
      <value value="89.1021788475929"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-ler">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-fert-cutoff">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-col">
      <value value="9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Interference_Sim_means" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1441"/>
    <metric>r3C</metric>
    <metric>r6C</metric>
    <metric>r9C</metric>
    <metric>r24C</metric>
    <metric>r3L</metric>
    <metric>r6L</metric>
    <metric>r9L</metric>
    <metric>r24L</metric>
    <metric>count patches with [pcolor = yellow]</metric>
    <metric>count patches with [pcolor = blue] / 6</metric>
    <metric>count patches with [pcolor = red] / 6</metric>
    <enumeratedValueSet variable="step-size">
      <value value="0.66"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-ler">
      <value value="18.975313187193"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-col">
      <value value="3.46683433710825"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accession">
      <value value="&quot;both&quot;"/>
      <value value="&quot;ler-only&quot;"/>
      <value value="&quot;col-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-type">
      <value value="&quot;best_mean_params&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-attrac-sat-fcn">
      <value value="&quot;Holling-Type-1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-starch">
      <value value="92.6892106482177"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-col">
      <value value="2.82183789445632"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-col">
      <value value="0.757218941401608"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-col">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="H1P1">
      <value value="5.34596861598446"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-ler">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-ler">
      <value value="5.64709511672291"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="LHS-paramcombo-counter" first="0" step="1" last="99"/>
    <enumeratedValueSet variable="prop-ler">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replenish-style">
      <value value="&quot;replenish-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-ler">
      <value value="89.1021788475929"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-ler">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-fert-cutoff">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-col">
      <value value="9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Interference_Sim_distr" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1441"/>
    <metric>r3C</metric>
    <metric>r6C</metric>
    <metric>r9C</metric>
    <metric>r24C</metric>
    <metric>r3L</metric>
    <metric>r6L</metric>
    <metric>r9L</metric>
    <metric>r24L</metric>
    <metric>count patches with [pcolor = yellow]</metric>
    <metric>count patches with [pcolor = blue] / 6</metric>
    <metric>count patches with [pcolor = red] / 6</metric>
    <enumeratedValueSet variable="step-size">
      <value value="0.66"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-ler">
      <value value="18.975313187193"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-col">
      <value value="3.46683433710825"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accession">
      <value value="&quot;both&quot;"/>
      <value value="&quot;ler-only&quot;"/>
      <value value="&quot;col-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-type">
      <value value="&quot;best_distr_params&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-attrac-sat-fcn">
      <value value="&quot;Holling-Type-1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-starch">
      <value value="92.6892106482177"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-col">
      <value value="2.82183789445632"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-col">
      <value value="0.757218941401608"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-col">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="H1P1">
      <value value="5.34596861598446"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-ler">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-ler">
      <value value="5.64709511672291"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="LHS-paramcombo-counter" first="0" step="1" last="99"/>
    <enumeratedValueSet variable="prop-ler">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replenish-style">
      <value value="&quot;replenish-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-ler">
      <value value="89.1021788475929"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-ler">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-fert-cutoff">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-col">
      <value value="9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="LerFertComp_distr" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1441"/>
    <metric>r3C</metric>
    <metric>r6C</metric>
    <metric>r9C</metric>
    <metric>r24C</metric>
    <metric>r3L</metric>
    <metric>r6L</metric>
    <metric>r9L</metric>
    <metric>r24L</metric>
    <metric>count patches with [pcolor = yellow]</metric>
    <metric>count patches with [pcolor = blue] / 6</metric>
    <metric>count patches with [pcolor = red] / 6</metric>
    <enumeratedValueSet variable="step-size">
      <value value="0.66"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-ler">
      <value value="18.975313187193"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-col">
      <value value="3.46683433710825"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accession">
      <value value="&quot;both&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-type">
      <value value="&quot;best_distr_params&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-attrac-sat-fcn">
      <value value="&quot;Holling-Type-1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-starch">
      <value value="92.6892106482177"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-col">
      <value value="2.82183789445632"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-col">
      <value value="0.757218941401608"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-col">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="H1P1">
      <value value="5.34596861598446"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-ler">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-ler">
      <value value="5.64709511672291"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="LHS-paramcombo-counter" first="0" step="1" last="99"/>
    <steppedValueSet variable="prop-ler" first="0.5" step="0.02" last="1"/>
    <enumeratedValueSet variable="replenish-style">
      <value value="&quot;replenish-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-ler">
      <value value="89.1021788475929"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-ler">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-fert-cutoff">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-col">
      <value value="9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="LerFertComp_means" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1441"/>
    <metric>r3C</metric>
    <metric>r6C</metric>
    <metric>r9C</metric>
    <metric>r24C</metric>
    <metric>r3L</metric>
    <metric>r6L</metric>
    <metric>r9L</metric>
    <metric>r24L</metric>
    <metric>count patches with [pcolor = yellow]</metric>
    <metric>count patches with [pcolor = blue] / 6</metric>
    <metric>count patches with [pcolor = red] / 6</metric>
    <enumeratedValueSet variable="step-size">
      <value value="0.66"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-ler">
      <value value="18.975313187193"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-col">
      <value value="3.46683433710825"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accession">
      <value value="&quot;both&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-type">
      <value value="&quot;best_mean_params&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-attrac-sat-fcn">
      <value value="&quot;Holling-Type-1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-starch">
      <value value="92.6892106482177"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-col">
      <value value="2.82183789445632"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert-thresh-col">
      <value value="0.757218941401608"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-col">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="H1P1">
      <value value="5.34596861598446"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-ler">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-starch-harvest-ler">
      <value value="5.64709511672291"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq-starch-replenish">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="LHS-paramcombo-counter" first="0" step="1" last="99"/>
    <steppedValueSet variable="prop-ler" first="0.5" step="0.02" last="1"/>
    <enumeratedValueSet variable="replenish-style">
      <value value="&quot;replenish-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-energy-ler">
      <value value="89.1021788475929"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chem-atr-rad-ler">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-fert-cutoff">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-cost-col">
      <value value="9"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
