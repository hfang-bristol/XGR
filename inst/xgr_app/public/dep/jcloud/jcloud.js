/**
 * jCloud
 * (C) 2011 Sawanna Team (http://sawanna.org)
 */


$.fn.jcloud=function(settings)
{
    jCloud={
        tags: new Array(),
        container_id: 'jcloud',
        radius: 100,
        size: 20,
        areal: 180,
        iterX: 0,
        iterY: 0,
        iterx: 1,
        itery: 1,
        step: 4,
        flats: 1,
        speed: 100,
        start_speed: 1,
        timer: 0,
        clock: 10,
        stop: true,
        splitX: 0,
        splitY: 0,
        colors: new Array(),
        init: function (elem,settings) {
            $(elem).css('display','none');
            
            var i=0;
            $(elem).find('li').each(function(){
                jCloud.tags[i]=$(this).html();
                i++;
            });
            
            if (typeof(settings) != 'undefined')
            {
                try {
                    if (typeof(settings.container_id) != 'undefined')
                    {
                        this.container_id=settings.container_id;
                    }
                    if (typeof(settings.radius) != 'undefined')
                    {
                        this.radius=settings.radius;
                    }
                    if (typeof(settings.size) != 'undefined')
                    {
                        this.size=settings.size;
                    }
                    if (typeof(settings.areal) != 'undefined')
                    {
                        this.areal=settings.areal;
                    }
                    if (typeof(settings.step) != 'undefined')
                    {
                        this.step=settings.step;
                    }
                    if (typeof(settings.flats) != 'undefined')
                    {
                        this.flats=settings.flats;
                    }
                    if (typeof(settings.speed) != 'undefined')
                    {
                        if (settings.speed < 1)
                        {
                            settings.speed=1;
                        }

                        this.speed=settings.speed;
                    }
                    if (typeof(settings.clock) != 'undefined')
                    {
                        this.clock=settings.clock;
                    }
                    if (typeof(settings.splitX) != 'undefined')
                    {
                        if (Math.abs(settings.splitX) < 2*this.radius)
                        {
                            this.splitX=settings.splitX;
                        }
                    }
                    if (typeof(settings.splitY) != 'undefined')
                    {
                        if (Math.abs(settings.splitY) < 2*this.radius)
                        {
                            this.splitY=settings.splitY;
                        }
                    }
                    if (typeof(settings.colors) != 'undefined' && settings.colors instanceof Array)
                    {
                        this.colors=settings.colors;
                    }
                } catch (e) {}
            }
            
            $(elem).replaceWith('<div id="'+this.container_id+'"></div>');
            
            $('#'+this.container_id).css({
                        'position':'relative',
                        'width':this.radius*2,
                        'height':this.radius*2,
                        'overflow':'hidden',
                        'padding': '0px'
                        });
                        
            $('#'+this.container_id).mousemove(function(e){
                var centerX=$(this).offset().left+jCloud.radius;
                var centerY=$(this).offset().top+jCloud.radius;
                
                jCloud.speed=jCloud.start_speed;
                jCloud.iterx=(centerX-e.pageX)/jCloud.radius;
                jCloud.itery=(centerY-e.pageY)/jCloud.radius;
            });
            
            $('#'+this.container_id).mouseover(function()
            {
                jCloud.stop=false;
            });
            
            $('#'+this.container_id).mouseout(function()
            {
                jCloud.stop=true;
            });
            
            $('body').append('<div id="jcloud-buffer"></div>');
            $('#jcloud-buffer').css({
                        'display': 'none',
                        'padding': '0px',
                        'white-space': 'nowrap'
                    });
        },
        render: function() {
            $('#'+this.container_id).empty();
            
            this.start_speed=this.speed;
            
            for (var j in this.tags)
            {
                this.draw(j);
            }
            
            this.cloud();
        },
        draw: function(id) {
            $('#jcloud-buffer').html('<div class="jtag" id="jtag-'+id+'">'+this.tags[id]+'</div>');
            
            var taghtml=$('#jcloud-buffer').html();
            var width=$('#jcloud-buffer').width();
            var height=$('#jcloud-buffer').height();
            
            $('#'+this.container_id).append(taghtml);
            
            $('#jtag-'+id).attr('width',width);
            $('#jtag-'+id).attr('height',height);
        },
        cloud: function(iter,flat_iter,update) {
            var count=this.tags.length;
            var per_flat=Math.ceil(count/this.flats);
            var delta=Math.round(360/per_flat);
            var deg=0;
            var flat_deg=0;

            if (typeof(iter) != 'undefined')
            {
                deg=Math.floor(iter);
            }
            
            if (typeof(flat_iter) != 'undefined')
            {
                flat_deg=Math.floor(flat_iter);
            }
            
            var flat_offset=Math.floor(this.areal/(this.flats+1))/2+flat_deg;
            var flat_delta=Math.floor(this.areal/this.flats);
            
            var flat_co=0;
            for (var j in this.tags)
            {
                this.pos(deg,j,parseInt($('#jtag-'+j).attr('width')),parseInt($('#jtag-'+j).attr('height')),flat_offset,update);
                deg+=delta;
                
                flat_co++;
                
                if (flat_co>=per_flat)
                {
                    flat_offset+=flat_delta;
                    flat_co=0;
                    deg+=flat_delta;
                }
            }
        },
        pos: function(degree,id,width,height,offset,update) {
            var radian=degree*Math.PI/180;
            var ell_offset=offset*Math.PI/180;
            var size_offset=0;
            
            var size=this.size+size_offset;
            
            var X=this.radius+Math.sin(radian)*this.radius*Math.cos(ell_offset)-width/2;
            var Y=this.radius-Math.cos(radian)*this.radius-height/2;
            
            size=(1+Math.cos(radian-1.57)*Math.sin(ell_offset))*size/2;
            
            if (size < 1)
            {
                size=1;
            }
            
            if (typeof(update) == 'undefined' || update==false)
            { 
                if (this.colors.length > 0)
                {
                    for (var i in this.colors)
                    {
                        if ((parseInt(id)+1) % (this.colors.length-i) == 0)
                        {
                            //$('#jtag-'+id).css('color',this.colors[this.colors.length-i-1]);
                            //$('#jtag-'+id).children('a').css('color',this.colors[this.colors.length-i-1]);
                            
                            $('#jtag-'+id).attr('style','color:'+this.colors[this.colors.length-i-1]+' !important');
                            $('#jtag-'+id).children('a').attr('style','color:'+this.colors[this.colors.length-i-1]+' !important');
                            break;
                        }
                    }
                }
                
                $('#jtag-'+id).css({
                    'position':'absolute',
                    'white-space': 'nowrap',
                    'font-weight': 'normal'
                });
            }
            
            var opacity=1;
            opacity=(1+Math.cos(radian-1.57)*Math.sin(ell_offset))*opacity/2;
            
            if (opacity<.5)
            {
                opacity=.5;
            }
            
            if (this.splitX)
            {
                var offsetX=parseInt(this.splitX);
                if (id%2==0)
                {
                    X-=offsetX;
                } else
                {
                    X+=offsetX;
                }
            } 
            
            if (this.splitY)
            {
                var offsetY=parseInt(this.splitY);
                if (id%2==0)
                {
                    Y-=offsetY;
                } else
                {
                    Y+=offsetY;
                }
            } 
            
            $('#jtag-'+id).css({
                    'left':Math.round(X),
                    'top':Math.round(Y),
                    'font-size':Math.floor(size)+'px',
                    'opacity': opacity
                });
            
            
                
        },
        update: function() {
            this.timer++;
            
            var interval=Math.floor(5/jCloud.speed);
            
            if (this.timer < interval)
            {
                return;
            }
            
            if (this.stop)
            {
                this.speed=this.speed-(this.speed*.1);
            }
            
            this.timer=0;
            
            if (this.iterX >= 360)
            {
                this.iterX=360-this.iterX;
            }
            
            if (this.iterY >= 360)
            {
                this.iterY=360-this.iterY;
            }
            
            this.iterX+=this.iterx;
            this.iterY+=this.itery;

            this.cloud(this.iterY,this.iterX,true);
        }
    }
    
    jCloud.init($(this),settings);
    
    jCloud.render();
    
    jCloudInterval=window.setInterval("jCloud.update()",jCloud.clock);
}
