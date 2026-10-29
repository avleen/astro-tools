<a href="../README.md"><img src="../docs/logo/logo.svg" width="64" align="right" alt="astro-tools"></a>

# nina-sequences

N.I.N.A. Advanced Sequencer files for the remote rig.

## Nightly - Target Scheduler.json

One sequence that runs unattended every night, for as long as N.I.N.A. is left
running. The [Target Scheduler](https://tcpalmer.github.io/nina-scheduler/)
plugin picks the targets; the sequence only handles the roof, start-up,
shut-down and flats. Built with N.I.N.A. 3.2 and Target Scheduler 5.9.

```
Nightly loop (repeats every night, up to 1000 times)
  1. Wait for roof open, then nautical dusk
       wait until safe · cool camera to -10 °C · wait for nautical dusk
  2. Night session, until nautical dawn or no targets left tonight
     2a. wait until safe · unpark · open cover · autofocus · start guiding
     2b. Target Scheduler, while safe, targets remain and before nautical dawn
     2c. stop guiding · find home · close cover
     (if the roof closes, 2b stops, 2c sends the mount home, and 2a waits for it to reopen)
  3. End of night
       stop guiding · find home · close cover
       Target Scheduler flats (25 per filter) · warm camera
  4. Wait until an hour after sunrise, then start the next night

Global triggers: meridian flip; autofocus after a 2 °C temperature change.
```

The flats in step 3 are what [masterflats](../masterflats) turns into master
flats each morning. See the [workflow example](../README.md#workflow-example).

### Using it

1. Install the Target Scheduler plugin and set up your projects and exposure
   plans in it. The targets aren't in this file.
2. In N.I.N.A., **Advanced Sequencer → Open** and load the file.
3. Check the camera temperature in step 1, the flat count in step 3 and the
   autofocus trigger, then start the sequence.

The rig needs a safety monitor (the roof), a cover or flat panel, and a guider.
The clock times shown in the waits (19:46, 06:56, ...) are just the last values
N.I.N.A. computed. Each wait follows its dusk, dawn or sunrise setting, so they
update for your site and date when the sequence loads.
