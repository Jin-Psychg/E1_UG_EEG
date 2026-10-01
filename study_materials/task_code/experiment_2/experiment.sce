############################################################
# Experiment 2 task scenario
# Identifying header comments removed for anonymous review.
# Task logic is otherwise unchanged from the archived source.
############################################################
 
no_logfile = true;  
active_buttons = 4;		
button_codes = 0,0,0,0; # 1,2 Tasten, 3 space, 4 mouse0
default_font_size = 24;
default_font = "arial";	
default_background_color = 170,143,127; # brown
default_text_color = 0,0,0; # black
default_formatted_text = true; # for HTML-tags (<,>,& sind Sonderzeichen)
response_matching = simple_matching;
#write_codes = true; # Trigger senden
pulse_width = 2;

####### SDL ##############################################################
begin;

picture {} default; # leerer Bildschirm

picture {box{height= 2; width=20; color=0,0,0;}; x=0; y=0; # ein Fixkreuz, nein - DAS Fixkreuz
			box{height=20; width= 2; color=0,0,0;}; x=0; y=0;}fk;
			
picture {}skala;

text {caption="......"; font_size=24; font="Arial"; max_text_width=1000;}frage; # für die scala
text {caption="links";  font_size=16; font="Arial"; max_text_width=450;}links; 
text {caption="rechts"; font_size=16; font="Arial"; max_text_width=450;}rechts; 
text {caption="mitte";  font_size=16; font="Arial";}mitte; 

picture {bitmap {filename="scale.png"; width=608;}verlauf; x=0; y=0;}verl; # Verlauf unter der ratigskala, verlauf.bmp

#************ Fixkreuze *****************************

# Fixkreuz
trial {trial_type = fixed;
		 trial_duration = 500;
		 picture fk;
		}fixkreuz;

#*********** Instruktion und Endebildschirm ****************************

trial {trial_duration=forever;
       trial_type=first_response;
		 picture {bitmap {preload=false;}instr; x=0; y=0;}pic1;			
      }instruk;		

trial {trial_duration=forever;
       trial_type=specific_response;
		 terminator_button = 1,2;      
		 picture {bitmap {preload=false;}instr_leader; x=0; y=0;}pic2;			
      }instruk_leader;	

trial {trial_type = fixed;
       trial_duration = 1500;
       picture {text{caption = "Ende des Hauptexperimentes";}; x=0; y=0;};
		}endescreen;
	
trial {trial_duration=forever;
       trial_type=first_response;	
		 picture {text{caption = "Pause\n- weiter mit Tastendruck -";}; x=0; y=0;};
		}pause;
		
#*********** Hauptroutinen ********************************************

trial {trial_type=first_response;
		 trial_duration = 3000;
	  	 stimulus_event{
			picture {bitmap {preload=false; width=200; scale_factor=scale_to_width;}pic; x=0; y=0;
						text{caption="o"; preload=false; font_size=14; font_color=255,255,0; trans_src_color=170,143,127;}oben;  x=0; y= 8;
						text{caption="u"; preload=false; font_size=14; font_color=255,255,0; trans_src_color=170,143,127;}unten; x=0; y=-7;
						text{caption="ar links";  preload=false; font_size=16;}ar_l; x=-62; y=-110;
						text{caption="ar rechts"; preload=false; font_size=16;}ar_r; x= 62; y=-110;};
			response_active=true;}pictrigger;
		}bild;
		
trial {trial_type = fixed;
		 trial_duration = 1000;
	  	 stimulus_event{
			picture {text{caption="o"; preload=false; font_size=20;}oben2;  x=0; y= 20;
						text{caption="u"; preload=false; font_size=20;}unten2; x=0; y=-20;
			}fb2;
		 response_active=true;}feedbacktrigger2;			
		}feedback2;

trial {trial_type = fixed;
		 trial_duration = 1000;
		  	 stimulus_event{
			picture {text{caption="zu langsam"; font_size=20;}; x=0; y=0;};
			port_code=100;}missedtrigger;
		}missed;		
	
###  PCL ############################################################
begin_pcl;

mouse mouse1 = response_manager.get_mouse(1);

include "arrays.pcl";

###Name des Logfiles eingeben
preset int subject = 0;
preset int blockstart = 0;
output_file file = new output_file;	 	

if blockstart!=0 then
 file.open(string(subject) + "_from" + string(blockstart) + ".txt", false);
else
 if subject == 0 then file.open(string(subject) + ".txt");
 else file.open(string(subject) + ".txt", false) end;
end;

if subject==0 then subject=1 end;

array <string> stim[0][0]; array <string> tmp[0][0];
string setting;

if 	 subject%10==1 then tmp.assign(order1); setting="1_A";
elseif subject%10==2 then tmp.assign(order2); setting="1_B";
elseif subject%10==3 then tmp.assign(order3); setting="2_A";
elseif subject%10==4 then tmp.assign(order4); setting="2_B";
elseif subject%10==5 then tmp.assign(order5); setting="3_A";
elseif subject%10==6 then tmp.assign(order6); setting="3_B";
elseif subject%10==7 then tmp.assign(order7); setting="4_A";
elseif subject%10==8 then tmp.assign(order8); setting="4_B";
elseif subject%10==9 then tmp.assign(order9); setting="5_A";
elseif subject%10==0 then tmp.assign(order10);setting="5_B" end;

int i;
int key=0;
int RT=0;
int AcRe=0; # accept/reject=1, reject/accept=2
int YouOther=0; # you/other=1, other/you=2
int reaction=0; # keycode
int prac=0;
int Rating; # 1-100 
int block;

#Kopfzeile in Logfile schreiben: \t = TAB, \n = Zeilenumbruch
file.print("subject\tsetting\tblock\tindex\tstim\tOffers_Other\tOffers_You\tOffer_Trigger\treaction\tRT\n");

# Logfile
sub logfile
begin
	file.print( string(subject)	+"\t"+
					setting				+"\t"+	
					string(block)		+"\t"+
					string(i)			+"\t"+
					stim[i][1]	 		+"\t"+
					stim[i][2]	 		+"\t"+	 		
					stim[i][3]	 		+"\t"+	 		
					stim[i][4]	 		+"\t"+	 		
					string(reaction)	+"\t"+
					string(RT)        +"\n");
end;

# Logfile
sub logfile2
begin
	file.print( string(subject)	+"\t"+
					setting				+"\t"+	
					string(block)		+"\t"+
					string(i)			+"\t"+
						 		"\t"+
						 		"\t"+	 		
						 		"\t"+	 		
						 		"\t"+	 		
					string(Rating)		+"\t"+
					string(RT)        +"\n");
end;


### instruktion
sub instruktion (string inhalt)
begin
	instr.set_filename (inhalt); instr.load();
instruk.present();
	instr.unload(); default.present(); wait_interval(750);
end;

### instruktion, die nur der Experimentleiter abbrechen kann
sub instruktion_leader (string inhalt)
begin
	instr_leader.set_filename (inhalt); instr_leader.load();
instruk_leader.present();
	instr_leader.unload(); default.present(); wait_interval(750);
end;

sub main (int start, int ende)
begin
 
	key=0; RT=0;
	if start==1 || start==281 || start==561 then AcRe=1 else AcRe=2 end; ### Annehmen <-> Ablehung blockweise wechseln 
	if prac==1 then missedtrigger.set_port_code(0) else missedtrigger.set_port_code(100) end;
		
	loop i=start until i>ende
	begin

	pic.set_filename(stim[i][1]+".png"); pic.load();

	YouOther=random(1,2);
	if YouOther==1 then oben.set_caption ("Sie  "+stim[i][3]+"0"); unten.set_caption("Part."+stim[i][2]+"0")
						else oben.set_caption ("Part."+stim[i][2]+"0"); unten.set_caption("Sie  "+stim[i][3]+"0") end;
	
	if AcRe==1 then ar_l.set_caption("Annehmen"); ar_l.set_font_color(59,125,35);  ar_r.set_caption("Ablehnung"); ar_r.set_font_color(192,0,0)
				  else ar_r.set_caption("Annehmen"); ar_r.set_font_color(59,125,35);  ar_l.set_caption("Ablehnung"); ar_l.set_font_color(192,0,0) end; 

	oben.load(); unten.load(); ar_l.load(); ar_r.load();
	pictrigger.set_port_code(int(stim[i][4])); 
	if prac==1 then pictrigger.set_port_code(0) end;
		
fixkreuz.present();
 if prac==0 then 
	if AcRe==1 then response_manager.set_button_codes({1,2,0,0}) else response_manager.set_button_codes({2,1,0,0}) end; ### Tasten-Trigger	
 end;
bild.present();

   response_manager.set_button_codes({0,0,0,0}); ### Tasten-Trigger aus
	stimulus_data last = stimulus_manager.last_stimulus_data();
	key = last.button(); RT = last.reaction_time();
	
	if key==AcRe then reaction=1 else reaction=2 end; if key==0 then reaction=0 end; # timeout
	
	if 	 reaction==1 && YouOther==1 then oben2.set_caption ("Sie erhalten <b>"+stim[i][3]+"0</b> ct"); unten2.set_caption("Part. erhält <b>"+stim[i][2]+"0</b> ct");
	elseif reaction==1 && YouOther==2 then oben2.set_caption ("Part. erhält <b>"+stim[i][2]+"0</b> ct"); unten2.set_caption("Sie erhalten <b>"+stim[i][3]+"0</b> ct"); 
	elseif reaction==2 && YouOther==1 then oben2.set_caption ("Sie erhalten <b>0</b> ct"); unten2.set_caption("Part. erhält <b>0</b> ct");	
	elseif reaction==2 && YouOther==2 then oben2.set_caption ("Part. erhält <b>0</b> ct"); unten2.set_caption("Sie erhalten <b>0</b> ct");
	end;	
	oben2.load(); unten2.load();

	if 	 stim[i][2]=="5" then feedbacktrigger2.set_port_code(3)
	elseif stim[i][2]=="6" then feedbacktrigger2.set_port_code(4)
	elseif stim[i][2]=="7" then feedbacktrigger2.set_port_code(5)
	elseif stim[i][2]=="8" then feedbacktrigger2.set_port_code(6)
	elseif stim[i][2]=="9" then feedbacktrigger2.set_port_code(7)
	end;

	if reaction==2 then feedbacktrigger2.set_port_code(8) end;
	if prac==1 then feedbacktrigger2.set_port_code(0) end;

wait_interval(100); ### triggerbreak

if key==0 then missed.present() else feedback2.present() end;

if prac==0 then logfile() end;

	pic.unload(); oben.unload(); unten.unload(); ar_l.unload(); ar_r.unload(); oben2.unload(); unten2.unload();
	
	if i%35==0 && i!=ende then instruktion("Folie7.PNG") end;

	i=i+1;
	end;
end;



### Rating
sub rating(string satz)
begin

##### Skala bauen ######
	polygon_graphic slider = new 	polygon_graphic; # Slider
	slider.set_line_width(10.0); slider.set_line_color(255,0,0,255); slider.set_sides(3); slider.set_rotation(180.0); slider.redraw(); skala.add_part(slider,0,0);

#	line_graphic slider = new line_graphic; # Slider
#	slider.set_line_width(10.0); slider.set_line_color(255,0,0,255); slider.add_line(0.0, -15.0, 0.0, 15.0); slider.redraw(); skala.add_part(slider,0,0);

	line_graphic track = new line_graphic; # Skala
	track.set_line_width(10.0); track.add_line(-300.0, 0.0, 300.0, 0.0); track.add_line(0.0, -20.0, 0.0, 20.0); track.set_line_color(170,143,127,255); track.redraw(); skala.add_part(track,0,0); # hidden
	
	mitte.set_caption("50%\nNeutral"); links.set_caption("0%\nÜberhaupt nicht"); rechts.set_caption("100%\nVollständig");
	links.redraw();  skala.add_part(links, -292, -70); # Beschriftung links
	rechts.redraw(); skala.add_part(rechts, 308, -70); # Beschriftung rechts	
	mitte.redraw();  skala.add_part(mitte, 12, -70);   # Beschriftung mitte	
	
	frage.set_caption(satz); frage.redraw(); skala.add_part(frage, 0, 200); # Frage hinzufügen

   skala.add_part(verlauf,0,30); # Verlauf hinzuf?gen
#### Skala-Ende ############

	skala.set_part_on_top(1,true); # Slider in den Vordergrund

	mouse1.set_min_max(1, -300,300); # Mausbereich festlegen
	mouse1.set_restricted(1, true);
	mouse1.set_xy(0,0);

#########################

int starttime=clock.time();

### Show scale until button 4 pressed
loop int count = response_manager.total_response_count(4) until response_manager.total_response_count(4) > count
begin
 mouse1.poll(); #read the mouse
 skala.set_part_x(1,mouse1.x());   #position the slider
 skala.present(); #draw it
end;

	RT=clock.time()-starttime; 

	Rating=((mouse1.x()/6)+50); #return rating 1-101; 
	if Rating==0 then Rating=1 end; 
	
###############
default.present(); wait_interval(1000); # ITI
logfile2();
end;

###############################################################
#################### sequence #################################
###############################################################

if blockstart==0 then

instruktion("Folie1.PNG"); # first instruction
instruktion("Folie2.PNG");
instruktion("Folie3.PNG");

prac=1; stim.assign(practice);
main(1,12); # practice

instruktion("Folie5.PNG");

prac=0; stim.assign(tmp);
instruktion("Folie3.PNG"); block=1; main( 1, 140); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=2; main(141,280); instruktion("Folie8.PNG");
instruktion("Folie3.PNG"); block=3; main(281,420); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=4; main(421,560); instruktion("Folie8.PNG");
instruktion("Folie3.PNG"); block=5; main(561,700); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=6; main(701,840);

instruktion("Folie9.PNG");

block=7;
i=1; rating("Inwieweit glauben Sie, dass das Foto des Gegenübers Ihre Entscheidung über die Vorschläge beeinflusst hat?");
i=2; rating("Inwieweit glauben Sie, dass die Vorschläge, auf die Sie gerade reagiert haben, von echten Menschen/Individuen gemacht wurden?");
i=3; rating("Inwieweit glauben Sie, dass die standardisierten Fotos auf echten Menschen/Individuen basieren?");

file.close();
instruktion("Folie10.PNG");
### end without blockselection ##############################


 elseif blockstart==1 then

prac=0; stim.assign(tmp);
instruktion("Folie3.PNG"); block=1; main( 1, 140); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=2; main(141,280); instruktion("Folie8.PNG");
instruktion("Folie3.PNG"); block=3; main(281,420); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=4; main(421,560); instruktion("Folie8.PNG");
instruktion("Folie3.PNG"); block=5; main(561,700); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=6; main(701,840);

instruktion("Folie9.PNG");

block=7;
i=1; rating("Inwieweit glauben Sie, dass das Foto des Gegenübers Ihre Entscheidung über die Vorschläge beeinflusst hat?");
i=2; rating("Inwieweit glauben Sie, dass die Vorschläge, auf die Sie gerade reagiert haben, von echten Menschen/Individuen gemacht wurden?");
i=3; rating("Inwieweit glauben Sie, dass die standardisierten Fotos auf echten Menschen/Individuen basieren?");

file.close();
instruktion("Folie10.PNG");
	
 elseif blockstart==2 then
	
prac=0; stim.assign(tmp);
instruktion("Folie4.PNG"); block=2; main(141,280); instruktion("Folie8.PNG");
instruktion("Folie3.PNG"); block=3; main(281,420); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=4; main(421,560); instruktion("Folie8.PNG");
instruktion("Folie3.PNG"); block=5; main(561,700); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=6; main(701,840);

instruktion("Folie9.PNG");

block=7;
i=1; rating("Inwieweit glauben Sie, dass das Foto des Gegenübers Ihre Entscheidung über die Vorschläge beeinflusst hat?");
i=2; rating("Inwieweit glauben Sie, dass die Vorschläge, auf die Sie gerade reagiert haben, von echten Menschen/Individuen gemacht wurden?");
i=3; rating("Inwieweit glauben Sie, dass die standardisierten Fotos auf echten Menschen/Individuen basieren?");

file.close();
instruktion("Folie10.PNG");

 elseif blockstart==3 then

prac=0; stim.assign(tmp);
instruktion("Folie3.PNG"); block=3; main(281,420); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=4; main(421,560); instruktion("Folie8.PNG");
instruktion("Folie3.PNG"); block=5; main(561,700); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=6; main(701,840);

instruktion("Folie9.PNG");

block=7;
i=1; rating("Inwieweit glauben Sie, dass das Foto des Gegenübers Ihre Entscheidung über die Vorschläge beeinflusst hat?");
i=2; rating("Inwieweit glauben Sie, dass die Vorschläge, auf die Sie gerade reagiert haben, von echten Menschen/Individuen gemacht wurden?");
i=3; rating("Inwieweit glauben Sie, dass die standardisierten Fotos auf echten Menschen/Individuen basieren?");

file.close();
instruktion("Folie10.PNG");

 elseif blockstart==4 then

prac=0; stim.assign(tmp);
instruktion("Folie4.PNG"); block=4; main(421,560); instruktion("Folie8.PNG");
instruktion("Folie3.PNG"); block=5; main(561,700); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=6; main(701,840);

instruktion("Folie9.PNG");

block=7;
i=1; rating("Inwieweit glauben Sie, dass das Foto des Gegenübers Ihre Entscheidung über die Vorschläge beeinflusst hat?");
i=2; rating("Inwieweit glauben Sie, dass die Vorschläge, auf die Sie gerade reagiert haben, von echten Menschen/Individuen gemacht wurden?");
i=3; rating("Inwieweit glauben Sie, dass die standardisierten Fotos auf echten Menschen/Individuen basieren?");

file.close();
instruktion("Folie10.PNG");

 elseif blockstart==5 then

prac=0; stim.assign(tmp);
instruktion("Folie3.PNG"); block=5; main(561,700); instruktion("Folie8.PNG");
instruktion("Folie4.PNG"); block=6; main(701,840);

instruktion("Folie9.PNG");

block=7;
i=1; rating("Inwieweit glauben Sie, dass das Foto des Gegenübers Ihre Entscheidung über die Vorschläge beeinflusst hat?");
i=2; rating("Inwieweit glauben Sie, dass die Vorschläge, auf die Sie gerade reagiert haben, von echten Menschen/Individuen gemacht wurden?");
i=3; rating("Inwieweit glauben Sie, dass die standardisierten Fotos auf echten Menschen/Individuen basieren?");

file.close(); 
instruktion("Folie10.PNG");

 elseif blockstart==6 then

prac=0; stim.assign(tmp);
instruktion("Folie4.PNG"); block=6; main(701,840);

instruktion("Folie9.PNG");

block=7;
i=1; rating("Inwieweit glauben Sie, dass das Foto des Gegenübers Ihre Entscheidung über die Vorschläge beeinflusst hat?");
i=2; rating("Inwieweit glauben Sie, dass die Vorschläge, auf die Sie gerade reagiert haben, von echten Menschen/Individuen gemacht wurden?");
i=3; rating("Inwieweit glauben Sie, dass die standardisierten Fotos auf echten Menschen/Individuen basieren?");

file.close(); 
instruktion("Folie10.PNG");
	
end;
